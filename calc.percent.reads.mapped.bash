#!/bin/bash


function usage { 
	echo "
	Usage: `basename $0` -b input_bam_dir [-o ./stats.mapping.tab]

	Given a directory of BAM files, reports mapping statistics
	to a specified tab-delimited output file. Requires bamtools.
	"
	}

nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case $1 in
		-b | --bams)  # dir of BAM files
			BAMDIR=$2
			echo "    calculating mapping statistics on BAMs in: $2"
			shift 2
			;;
		-o | --out)  # output file (optionally with filepath)
			OUT=$2
			echo "    output file saved as: $2"
			shift 2
			;;
		-h | --help)
			usage
			exit 1
			;;
	esac
done

# depend checks
[[ -z $BAMDIR ]] && { usage; exit 1; }
command -v bamtools >/dev/null 2>&1 || { echo 'ERROR: bamtools not found' >&2; exit 1; }

# default output filename
[[ -z $OUT ]] && OUT='stats.mapping.tab'

# Create output path and output summary file
DIR=$(dirname "$OUT")
if [ ! -d "$DIR" ]; then
	mkdir -p "$DIR"
	echo '    Created output dir...'
else
	echo "WARNING: Output dir $DIR already exists. Competing files will be overwritten."
fi
echo -e "Sample\tReads_Mapped[%]" > "$OUT" # overwrite to avoid appending data to different projects (if script exec >1 with same outfile)

for bam in "$BAMDIR"/*.bam; do
	b=$(basename "$bam")
	bamtools stats -in "$bam" > "$DIR"/tmpstatsfull."$b".tmp
	
	grep 'Mapped reads' tmpstatsfull."$b".tmp |\
	 awk 'BEGIN {FS="\t"}; {print $2}' |\
	 sed 's/[()]//g' > "$DIR"/tmpstatsperc."$b".tmp
	rm "$DIR"/tmpstatsfull."$b".tmp

	if [ -s tmpstatsperc."$b".tmp ]; then # percent was extracted
		p=$(cat "$DIR"/tmpstatsperc."$b".tmp)
		echo -e "$b\t$p" >> "$OUT"
	else # no percent avail
		echo -e "$b\tError" >> "$OUT"
	fi
	rm "$DIR"/tmpstatsperc."$b".tmp
done

