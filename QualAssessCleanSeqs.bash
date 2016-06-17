#!/bin/bash
#  Dependencies:  awk, bc, cut, GNU sed, grep, gunzip, mkdir, readlink, paste, printf, wc

usage() { 
	echo "
	USAGE: $0 Input_R1.paired.fq Input_R2.paired.fq Input.single.fq Outdir

	The trimmed reads must be in Illumina 1.8+ FastQ
	format, and filenames must be formatted:
	<name>.fastq.gz or <name>.fastq, where
	<name> is the name of the sample.
	"
	}

#Help and usage information
if [[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]]; then
	usage
	exit 1
elif [[ "$1" == "-v" || "$1" == "--version" ]]; then
	echo 'v1.1'
	exit 0
fi

#Require 4 arguments
if [ $# -ne 4 ]; then
	usage
	exit 1
fi

#Make outdir if absent
if [ ! -d "$4" ]; then
	mkdir -p "$4"
	echo "created output directory path:  $4"
fi

#File input sequence information
#  file regex from Illumina 1.8+ (gunzip extracted):
#  ([0-9A-Za-z]{1,19})_S([0-9A-Za-z]{1,3})_L([0-9]{3})_R[12]_([0-9]{3}).fastq
SAMPLE=$(basename $1 | cut -d . -f 1 | cut -d _ -f 1)

#Filename handling
if [[ $1 == *.gz && $2 == *.gz && $3 == *.gz ]]; then  #all files are gunzipped
	SAMPLEL=$(basename $1 .gz)
	SAMPLER=$(basename $2 .gz)
	SAMPLES=$(basename $3 .gz)
	echo "extracting the gunzipped file $1..."
	gunzip -c "$1" > "$4"/"$SAMPLEL"
	echo "$1 was extracted..."
	SEQ1="$4"/"$SAMPLEL"
	echo "extracting the gunzipped file $2..."
	gunzip -c "$2" > "$4"/"$SAMPLER"
	echo "$2 was extracted..."
	SEQ2="$4"/"$SAMPLER"
	echo "extracting the gunzipped file $3..."
	gunzip -c "$3" > "$4"/"$SAMPLES"
	echo "$3 was extracted..."
	SEQ3="$4"/"$SAMPLES"
elif [[ $1 == *q && $2 == *q && $3 == *q ]]; then  #all files are fq format
	echo "all three files are in FastQ format..."
	SEQ1="$1"
	SEQ2="$2"
	SEQ3="$3"
else  #for now skip possibily of some but not all gunzipped
	usage
	exit 1
fi
wait

echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " $SEQ1, $SEQ2,"
echo " and $SEQ3 sequences were read in..."
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Count total nucleotides
TN1=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ1)
TN2=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ2)
TN3=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ3)
BP_TOT="$(($TN1 + $TN2 + $TN3))"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " counted the total nucleotides..."
echo "        Total:  $BP_TOT"
echo "        R1:     $TN1"
echo "        R2:     $TN2"
echo "        single: $TN3"
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

OPSYS=$(uname -s)

###Count nucleotides >=Q20 and >=Q30 for Illumina 1.8+
if [[ "$OPSYS" == Linux ]]; then  #should be linux where sed default is GNU sed
	Q301=$(sed -n '4~4'p $SEQ1 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q302=$(sed -n '4~4'p $SEQ2 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q303=$(sed -n '4~4'p $SEQ3 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q20to301=$(sed -n '4~4'p $SEQ1 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m)
	Q20to302=$(sed -n '4~4'p $SEQ2 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m) 
	Q20to303=$(sed -n '4~4'p $SEQ3 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m) 
	Q201=$(($Q20to301 + $Q301))
	Q202=$(($Q20to302 + $Q302)) 
	Q203=$(($Q20to303 + $Q303)) 
elif [[ "$OPSYS" == Darwin ]]; then  #use gsed instead of sed in Mac OS X here
	Q301=$(gsed -n '4~4'p $SEQ1 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q302=$(gsed -n '4~4'p $SEQ2 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q303=$(gsed -n '4~4'p $SEQ3 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q20to301=$(gsed -n '4~4'p $SEQ1 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m)
	Q20to302=$(gsed -n '4~4'p $SEQ2 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m) 
	Q20to303=$(gsed -n '4~4'p $SEQ3 | grep -o '[5-9:\;\<=\>]' | paste -s -d"\0" - | wc -m) 
	Q201=$(($Q20to301 + $Q301))
	Q202=$(($Q20to302 + $Q302)) 
	Q203=$(($Q20to303 + $Q303)) 
else
	echo 'ERROR: Your operating system does not appear to be supported.' >&2
	exit 1
fi
wait
Q20_TOT="$(($Q201 + $Q202 + $Q203))" 
Q30_TOT="$(($Q301 + $Q302 + $Q303))" 
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' counted the Q20+ and Q30+ nucleotides...'
echo "        ${Q201} nucleotides are >=Q20 in $SEQ1" 
echo "        ${Q301} nucleotides are >=Q30 in $SEQ1" 
echo "        ${Q202} nucleotides are >=Q20 in $SEQ2" 
echo "        ${Q302} nucleotides are >=Q30 in $SEQ2" 
echo "        ${Q203} nucleotides are >=Q20 in $SEQ3" 
echo "        ${Q303} nucleotides are >=Q30 in $SEQ3" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Calculate percentages of >=Q20 and >=Q30
PQ201="$(echo "scale=2;($Q201/$TN1)*100" | bc)"
PQ202="$(echo "scale=2;($Q202/$TN2)*100" | bc)"
PQ203="$(echo "scale=2;($Q203/$TN3)*100" | bc)"
PQ301="$(echo "scale=2;($Q301/$TN1)*100" | bc)"
PQ302="$(echo "scale=2;($Q302/$TN2)*100" | bc)"
PQ303="$(echo "scale=2;($Q303/$TN3)*100" | bc)"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' calculated the percentages of Q20 and Q30...'
echo "        ${PQ201}% nucleotides are >=Q20 in $SEQ1" 
echo "        ${PQ301}% nucleotides are >=Q30 in $SEQ1" 
echo "        ${PQ202}% nucleotides are >=Q20 in $SEQ2" 
echo "        ${PQ302}% nucleotides are >=Q30 in $SEQ2" 
echo "        ${PQ203}% nucleotides are >=Q20 in $SEQ3" 
echo "        ${PQ303}% nucleotides are >=Q30 in $SEQ3" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Count number of reads
NUMLN_L=$(wc -l $SEQ1| awk '{print $1}')
NUMLN_R=$(wc -l $SEQ2 | awk '{print $1}')
NUMLN_single=$(wc -l $SEQ3 | awk '{print $1}')
wait
NUM_L="$(($NUMLN_L / 4))"
NUM_R="$(($NUMLN_R / 4))"
NUM_single="$(($NUMLN_single / 4))"
READ_TOT="$(($NUM_L + NUM_R + NUM_single))"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' counted the number of reads...'
echo "        There are $NUM_L reads in $SAMPLEL" 
echo "        There are $NUM_R reads in $SAMPLER" 
echo "        There are $NUM_single reads in $SAMPLES" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

#Create results output file
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
Sample_Name Q20_Total_[bp] Q30_Total_[bp] Q20_R1_[bp] Q20_R2_[bp] Q20_single_[bp] Q20_R1_[%] Q20_R2_[%] Q20_single_[%] Q30_R1_[bp] Q30_R2_[bp] Q30_single_[bp] Q30_R1_[%] Q30_R2_[%] Q30_single_[%] Total_Sequenced_[bp] Total_Sequenced_[reads] NUMR1_[bp] NUMR2_[bp] NUMsingle_[bp] \
"$SAMPLE" "$Q20_TOT" "$Q30_TOT" "$Q201" "$Q202" "$Q203" "$PQ201%" "$PQ202%" "$PQ203%" "$Q301" "$Q302" "$Q303" "$PQ301%" "$PQ302%" "$PQ303%" "$BP_TOT" "$READ_TOT" "$NUM_L" "$NUM_R" "$NUM_single" 
> "$4"/QualAssessTrimSeqs_"$SAMPLE"_results.tab ;
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' created results output file...'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

if [[ "$OPSYS" == Linux ]]; then
	R1FILEPATH=$(readlink -m $1)
	R2FILEPATH=$(readlink -m $2)
	singleFILEPATH=$(readlink -m $3)
elif [[ "$OPSYS" == Darwin ]]; then
	R1FILEPATH=$(echo `pwd`/`ls "$1"`)
	R2FILEPATH=$(echo `pwd`/`ls "$2"`)
	singleFILEPATH=$(echo `pwd`/`ls "$3"`)
fi

[[ $1 == *.gz ]] && rm -f "$3"/"$SAMPLEL" 
[[ $2 == *.gz ]] && rm -f "$3"/"$SAMPLER"
[[ $3 == *.gz ]] && rm -f "$3"/"$SAMPLES"

#Create log file
printf "`date`\n$USER\n%s\n%s\n\n" \
"$R1FILEPATH" "$R2FILEPATH" "$singleFILEPATH" \
> "$4"/QualAssessTrimSeqs_"$SAMPLE"_results.log ;
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' created log output file...'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " Quality assessment of $SAMPLE trimmed sequences completed" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
