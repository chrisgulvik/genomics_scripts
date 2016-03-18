#!/bin/bash

set -e
function usage() { 
	echo "
	Given a FastA file, identifies repetitive regions, and outputs a BED file

	usage: `basename $0` -r reference.fasta -o maskedRegions.bed [-l <int>] [-i <float>] [-t <dir>]
		
		-r [default: ./reference.fasta]   FastA formatted input file
		-o [default: ./maskedRegions.bed] BED formatted output file
		-i [default: 99]                  Minimum percent identity required between two
		                                  repeat regions
		-l [default: 249]                 Minimum repeat length to find (in bp)
		-t [default: ./tmp/]              Temporary directory where intermediate files
		                                  are written to and removed
		
	Paths can be absolute or relative.
	"
	}

# Defaults
PERC_ID=99
LEN=249
OUTFILE="$PWD/maskedRegions.bed"
REF="$PWD/reference.fasta"
TMP="$PWD/tmp"

nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case "$1" in
		-i | --percid)
			PERC_ID="$2"
			echo "    min percent identity required: $2"
			shift 2
			;;
		-l | --len)
			LEN="$2"
			echo "    min length: $2"
			shift 2
			;;
		-o | --out)
			OUTFILE="$2"
			echo "    output file: $2"
			shift 2
			;;
		-r | --ref)
			REF="$2"
			echo "    reference genome: $2"
			shift 2
			;;
		-t | --temp)
			TMP="$2"
			echo "    temporary directory: $2"
			shift 2
			;;
		-h | --help)
			usage
			exit 1
			;;
		\?)
			echo "    ERROR: $2 is not a valid argument"
			usage
			exit 1
			;;
	esac
done

# Input file requirements and dependency check
[[ ! -e "$REF" || ! -s "$REF" ]] && { echo "ERROR: $REF cannot be read"; exit 1; }
command -v nucmer >/dev/null 2>&1 || { echo 'ERROR: nucmer not found' >&2; exit 1; }
command -v show-coords >/dev/null 2>&1 || { echo 'ERROR: show-coords not found' >&2; exit 1; }

# Create tmp dir (if absent)
if [ ! -d "$TMP" ]; then
	mkdir -p "$TMP"
else
	echo "WARNING: Temporary directory $TMP already exists."
	echo 'Files present in the specified tmp dir might disrupt analysis.'
fi

# Primary system commands
cd "$TMP"  #cannot specific outpath in nucmer, so change into tmp dir
nucmer --nosimplify --maxmatch -p DUPES "$REF" "$REF" 2> /dev/null
show-coords -I "$PERC_ID" -L "$LEN" -r -l -o -T DUPES.delta > filtered.coords

# Parse nucmer output into BED format
[[ -a DUPES.bed ]] && rm DUPES.bed
while IFS=$'\t' read -r -a line; do
	if [ ${line[0]} -gt ${line[1]} ]; then
		l0=${line[0]}
		l1=${line[1]}
		line[0]=$l1
		line[1]=$l0
	fi
	if [ ${line[2]} -gt ${line[3]} ]; then
		l2=${line[2]}
		l3=${line[3]}
		line[2]=$l3
		line[3]=$l2
	fi
	echo -e "${line[9]}\t${line[0]}\t${line[1]}" >> DUPES.bed
	echo -e "${line[9]}\t${line[2]}\t${line[3]}" >> DUPES.bed
done < <(grep '^[1-9]' filtered.coords)

sort -t$'\t' -k2 -n DUPES.bed > DUPES_sorted.bed
awk '{if(++dup[$0]==1) {print $0}}' DUPES_sorted.bed > DUPES_sorted_deduped.bed

# Remove full-length chromosome matches; handles >1 contig
while read -r -a deflines; do
	CHROM_ID=${deflines[1]}
	END_POS=${deflines[3]}
	# awk -v chrom="${CHROM_ID}" -v end="${END_POS}" '!/chrom\t1\tend/' DUPES_sorted_deduped.bed > tmp.bed
	# mv tmp.bed DUPES_sorted_deduped.bed
	TAB=$'\t'
	sed -i "/${CHROM_ID}${TAB}1${TAB}${END_POS}/d" DUPES_sorted_deduped.bed
done < <(grep '^>' DUPES.delta)

# Merge common overlapping regions
bedtools merge -i DUPES_sorted_deduped.bed > maskedRegions.bed

REGIONS=$(wc -l maskedRegions.bed | awk '{print $1}')
SITES=$(awk -F'\t' 'BEGIN{SUM=0}; {SUM+=$3-$2}; END{print SUM}' maskedRegions.bed)
echo "found $REGIONS regions and $SITES total sites to mask"

# Cleanup
mv -v maskedRegions.bed "$OUTFILE"
rm -r "$TMP"
