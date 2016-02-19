#!/bin/bash
#  Dependencies:  awk, bc, cut, GNU sed, grep, gunzip, mkdir, readlink, paste, printf, wc

usage() { 
	echo "
	USAGE: $0 Input_R1.fastq Input_R2.fastq Outdir

	The raw reads must be in Illumina 1.8+ FastQ
	format, and filenames must be formatted:
	<name>.fastq.gz or <name>.fastq, where
	<name> is the name of the sample.
	"
	}

#Help and usage information
if [[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]]; then
	usage
	exit 1
fi

#Require 3 arguments
if [ $# -ne 3 ]; then
	usage
	exit 1
fi

#Make outdir if absent
if [ ! -d "$3" ]; then
	mkdir -p "$3"
	echo "created output directory path:  $3"
fi

#File input sequence information
#  file regex from Illumina 1.8+ (gunzip extracted):
#  ([0-9A-Za-z]{1,19})_S([0-9A-Za-z]{1,3})_L([0-9]{3})_R[12]_([0-9]{3}).fastq
SAMPLE=$(basename $1 | cut -d . -f 1 | cut -d _ -f 1,2)
SAMPLEL=$(basename $1)
SAMPLER=$(basename $2)

#Filename handling
dIFS=$IFS  #default IFS
IFS=$''  #change for newlines
if [[ $1 == *.gz && $2 == *.gz ]]; then  #both files are gunzipped
	SAMPLEL=$(basename $SEQ1 .gz)
	SAMPLER=$(basename $SEQ2 .gz)
	echo "extracting the gunzipped file $1..."
	gunzip -c $1 > "$3"/"$SAMPLEL"
	echo "$1 was extracted..."
	SEQ1=$("$3"/"$SAMPLEL")
	echo "extracting the gunzipped file $2..."
	gunzip -c $2 > "$3"/"$SAMPLER"
	echo "$2 was extracted..."
	SEQ2=$("$3"/"$SAMPLER")
elif [[ $1 == *q && $2 == *q ]]; then  #both files are fastq format
	echo "both are FastQ format..."
	SEQ1=$(echo $1)
	SEQ2=$(echo $2)
elif [[ $1 == *.gz && $2 == *q ]]; then
	SAMPLEL=$(basename $SEQ1 .gz)
	echo "extracting the gunzipped file $1..."
	gunzip -c $1 > "$3"/"$SAMPLEL"
	echo "$1 was extracted..."
	SEQ1=$("$3"/"$SAMPLEL")
	SEQ2="$2"  #$2 is still in fastq format
elif [[ $1 == *q && $2 == *.gz ]]; then
	SAMPLER=$(basename $SEQ2 .gz)
	SEQ1="$1"  #$1 is still in fastq format
	echo "extracting the gunzipped file $2..."
	gunzip -c $2 > "$3"/"$SAMPLER"
	echo "$2 was extracted..."
	SEQ2=$("$3"/"$SAMPLER")
else
	usage
	exit 1
fi
IFS=$dIFS  #back to default IFS
wait

echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " $SAMPLEL and $SAMPLER sequences were read in..."
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Count total nucleotides
TN1=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ1)
TN2=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ2)
BP_TOT="$(($TN1 + $TN2))"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " counted the total nucleotides..."
echo "        Total: $BP_TOT"
echo "        R1:    $TN1"
echo "        R2:    $TN2"
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

OPSYS=$(uname -s)

#Count nucleotides >=Q20 and >=Q30 for Illumina 1.8+
if [[ "$OPSYS" == Linux ]]; then  #should be linux where sed default is GNU sed
	Q30L=$(sed -n '4~4'p $SEQ1 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q30R=$(sed -n '4~4'p $SEQ2 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q20to30L=$(sed -n '4~4'p $SEQ1 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m)
	Q20to30R=$(sed -n '4~4'p $SEQ2 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m) 
	Q20L="$(($Q20to30L + $Q30L))"
	Q20R="$(($Q20to30R + $Q30R))"
elif [[ "$OPSYS" == Darwin ]]; then  #use gsed instead of sed in Mac OS X here
	Q30L=$(gsed -n '4~4'p $SEQ1 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q30R=$(gsed -n '4~4'p $SEQ2 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q20to30L=$(gsed -n '4~4'p $SEQ1 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m)
	Q20to30R=$(gsed -n '4~4'p $SEQ2 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m) 
	Q20L="$(($Q20to30L + $Q30L))"
	Q20R="$(($Q20to30R + $Q30R))"
else
	echo 'ERROR: Your operating system does not appear to be supported.' >&2
	exit 1
fi
wait
Q20_TOT="$(($Q20L + $Q20R))" 
Q30_TOT="$(($Q30L + $Q30R))" 
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' counted the Q20+ and Q30+ nucleotides...'
echo "        ${Q20L} nucleotides are >=Q20 in $SAMPLEL" 
echo "        ${Q30L} nucleotides are >=Q30 in $SAMPLEL" 
echo "        ${Q20R} nucleotides are >=Q20 in $SAMPLER" 
echo "        ${Q30R} nucleotides are >=Q30 in $SAMPLER" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Calculate percentages of >=Q20 and >=Q30
PQ20L="$(echo "scale=2;($Q20L/$TN1)*100" | bc)"
PQ20R="$(echo "scale=2;($Q20R/$TN2)*100" | bc)"
PQ30L="$(echo "scale=2;($Q30L/$TN1)*100" | bc)"
PQ30R="$(echo "scale=2;($Q30R/$TN2)*100" | bc)"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' calculated the percentages of Q20+ and Q30+...'
echo "        ${PQ20L}% nucleotides are >=Q20 in $SAMPLEL" 
echo "        ${PQ30L}% nucleotides are >=Q30 in $SAMPLEL" 
echo "        ${PQ20R}% nucleotides are >=Q20 in $SAMPLER" 
echo "        ${PQ30R}% nucleotides are >=Q30 in $SAMPLER" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Count number of reads
NUMLN_L=$(wc -l $SEQ1| awk '{print $1}')
NUMLN_R=$(wc -l $SEQ2 | awk '{print $1}')
wait
NUM_L="$(($NUMLN_L / 4))"
NUM_R="$(($NUMLN_R / 4))"
READ_TOT="$(($NUM_L + NUM_R))"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' counted the number of reads...'
echo "        There are $NUM_L reads in $SAMPLEL" 
echo "        There are $NUM_R reads in $SAMPLER" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Create results output file
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
Sample_Name Q20_Total_[bp] Q30_Total_[bp] Q20_R1_[bp] Q20_R2_[bp] Q20_R1_[%] Q20_R2_[%] Q30_R1_[bp] Q30_R2_[bp] Q30_R1_[%] Q30_R2_[%] Total_Sequenced_[bp] Total_Sequenced_[reads] NUMR1_[bp] NUMR2_[bp] \
"$SAMPLE" "$Q20_TOT" "$Q30_TOT" "$Q20L" "$Q20R" "$PQ20L%" "$PQ20R%" "$Q30L" "$Q30R" "$PQ30L%" "$PQ30R%" "$BP_TOT" "$READ_TOT" "$NUM_L" "$NUM_R" > "$3"/QualAssessRawSeqs_"$SAMPLE"_results.tsv ;
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' created results output file...'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

if [[ "$OPSYS" == Linux ]]; then
	R1FILEPATH=$(readlink -m $1)
	R2FILEPATH=$(readlink -m $2)
elif [[ "$OPSYS" == Darwin ]]; then
	R1FILEPATH=$(echo `pwd`/`ls "$1"`)
	R2FILEPATH=$(echo `pwd`/`ls "$2"`)
fi

#Create log file
printf "`date`\n$USER\n%s\n%s\n\n" \
"$R1FILEPATH" "$R2FILEPATH" \
> "$3"/QualAssessRawSeqs_"$SAMPLE"_results.log ;
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' created log output file...'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " Quality assessment of $SAMPLE raw sequences completed" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
