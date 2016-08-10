#!/bin/bash

function usage { 
	echo "
	Usage: `basename $0` /InputPath/input.fasta
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
TOT_NUM_A=$(grep -o 'A' /tmp/sequence.txt | wc -l)
TOT_NUM_T=$(grep -o 'T' /tmp/sequence.txt | wc -l)
TOT_NUM_C=$(grep -o 'C' /tmp/sequence.txt | wc -l)
TOT_NUM_G=$(grep -o 'G' /tmp/sequence.txt | wc -l)
TOT_NUM_N=$(grep -o 'N' /tmp/sequence.txt | wc -l)
TOT_NUM_GAP=$(grep -o '-' /tmp/sequence.txt | wc -l)
rm /tmp/sequence.txt

echo "A: $TOT_NUM_A
T: $TOT_NUM_T
C: $TOT_NUM_C
G: $TOT_NUM_G
N: $TOT_NUM_N
gaps: $TOT_NUM_GAP"
