#!/bin/sh


[[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]] && {echo 'Proj dir required as input'; exit 1;}

# Setup Proj where RawFQs already exist
mkdir -p "$1"/{TrimFQs,Asm}
cd "$1"/RawFQs

# Qual Trim
for i in *_L001_R1_001.fastq; do
	b=$(basename $i _R1_001.fastq);
	SolexaQA++ dynamictrim $b_L001_R1_001.fastq $b_L001_R2_001.fastq -d ../TrimFQs --sanger -h 20;
	SolexaQA++ lengthsort ../TrimFQs/$b_L001_R1_001.fastq.trimmed ../TrimFQs/$b_L001_R2_001.fastq.trimmed -d ../TrimFQs -l 50;
done

# Merge PDF summaries of trimming
cd ../TrimFQs
gs -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -q -sOutputFile=Summary_Q20_Trim.pdf *.fastq_trimmed.segments_hist.pdf
gs -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -q -sOutputFile=Summary_Q20_Trim.Sort.pdf *.fastq.trimmed.summary.txt.pdf
rm *.fastq_trimmed.segments_hist.pdf *.fastq.trimmed.summary.txt.pdf 

# Rename file extensions because my QA script && SPAdes mandate FQ or FASTQ extension
rename "s/\.fastq\.trimmed\.paired/\.trimmed\.paired\.fq/" *.fastq.trimmed.paired
rename "s/\.fastq\.trimmed\.single/\.trimmed\.single\.fq/" *.fastq.trimmed.single

# Assemble
for i in *.trimmed.single.fq; do
	b=$(basename $i _L001_R1_001.trimmed.single.fq);
	spades.py -k 41,79,85,97 --careful --phred-offset 33 --pe1-1 "$b"_L001_R1_001.trimmed.paired.fq --pe1-2 "$b"_L001_R2_001.trimmed.paired.fq --pe1-s "$b"_L001_R1_001.trimmed.single.fq -o ../Asm/$b -t 30;
done
