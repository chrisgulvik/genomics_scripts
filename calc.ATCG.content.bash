#!/bin/bash


function usage { 
	echo "
	Usage: `basename $0` input.fasta

	Given a FastA file (with or without linewraps),
	reports nucleotide frequency including gaps and Ns,
	GC content, and cumulative length.
	"
	}

# check for infile
[[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]] && { usage; exit 1;}
if test -f "$1" -a -r "$1"; then
	GENOME=$1  # input FastA file
else
	usage; exit 1;
fi

# Calculate full ATCG content
grep -v '^>' "$GENOME" | tr -d [:space:] > /tmp/sequence.txt
TOT_NUM_A=$(grep -io 'A' /tmp/sequence.txt | wc -l)
TOT_NUM_T=$(grep -io 'T' /tmp/sequence.txt | wc -l)
TOT_NUM_C=$(grep -io 'C' /tmp/sequence.txt | wc -l)
TOT_NUM_G=$(grep -io 'G' /tmp/sequence.txt | wc -l)
TOT_NUM_N=$(grep -io 'N' /tmp/sequence.txt | wc -l)
TOT_NUM_GAP=$(grep -o '-' /tmp/sequence.txt | wc -l)
rm /tmp/sequence.txt

GC=$(($TOT_NUM_G + $TOT_NUM_C))
TOT_LEN=$(($GC + $TOT_NUM_A + $TOT_NUM_T + $TOT_NUM_N + $TOT_NUM_GAP))
GC_PERC=$(echo "scale=1; ($GC / $TOT_LEN) * 100" | bc)

echo "A: $TOT_NUM_A
T: $TOT_NUM_T
C: $TOT_NUM_C
G: $TOT_NUM_G
N: $TOT_NUM_N
gaps: $TOT_NUM_GAP
GC content: $GC_PERC%
Total length: $TOT_LEN
"
