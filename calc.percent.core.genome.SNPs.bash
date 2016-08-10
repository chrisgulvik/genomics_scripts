#!/bin/bash


function usage { 
	echo "
	Usage: `basename $0` -c coreSNPs.fasta -r reference.fasta
	"
	}

nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case $1 in
		-c | --core)  # FastA format
			CORESNPS=$2
			echo "    core SNP file: $2"
			shift 2
			;;
		-r | --ref)  # FastA format
			REF=$2
			echo "    ref genome file: $2"
			shift 2
			;;
		-h | --help)
			usage
			exit 1
			;;
	esac
done

[[ -z $CORESNPS || -z $REF ]] && { usage; exit 1; }

REF_GENOME=$(grep -v '^>' $REF | tr -d [:space:] | wc -c)
echo "    found $REF_GENOME sites in the reference genome"

TOT_SITES=$(grep -v '^>' $CORESNPS | tr -d [:space:] | wc -c)
TOT_SMPL=$(grep -c '^>' $CORESNPS)
SITES_PER_SMPL=$((TOT_SITES / TOT_SMPL))
echo "    found $SITES_PER_SMPL sites in each sample in $CORESNPS"

CORE_PER_REF=$(echo "scale=4;($SITES_PER_SMPL/$REF_GENOME)*100" | bc)
echo -e "\n    ${CORE_PER_REF}% of the reference genome was identified as core sites\n"
