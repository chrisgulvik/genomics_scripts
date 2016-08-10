#!/bin/bash

# recursively searches within input dir to remove most dirs and files SPAdes creates
# but maintains essential log files, contigs, and scaffolds FastA files

# usage: bash ~/scripts/prune_SPAdes_assembly_dirs.bash <PARENT_DIR_CLEANUP>

mapfile -t asm_dirs < <(find -L "$1" -type d -regextype posix-extended -regex '\/.*\/K[1-9][0-9]{1,2}' -print | xargs dirname | uniq)

for directory in "${asm_dirs[@]}"; do 
	rm -v "$directory"/{assembly_graph.fastg,before_rr.fasta,contigs.paths,dataset.info,input_dataset.yaml,scaffolds.paths};
done

find -L "$1" -type d -regextype posix-extended -regex '\/.*\/(K[1-9][0-9]{1,2}|misc|tmp)' -print | xargs rm -rv
