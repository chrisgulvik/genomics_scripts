#!/bin/bash


if [[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]]; then
  echo "
  Usage: `basename $0` input_dir

  Given a path containing SPAdes assemblies, recursively
  searches within the input dir to remove most dirs and files
  SPAdes creates but maintains essential contigs, graph, and
  log files."
  exit 1
fi

mapfile -t asm_dirs < <(find -L "$1" -type d \
 -regextype posix-extended -regex '\/.*\/K[1-9][0-9]{1,2}' -print |\
 xargs dirname | uniq)

for directory in "${asm_dirs[@]}"; do 
  rm -f "$directory"/{assembly_graph.fastg,before_rr.fasta}
  rm -f "$directory"/{contigs.paths,dataset.info,input_dataset.yaml}
  rm -f "$directory"/scaffolds.{fasta,paths}
  pigz -9f "$directory"/{params.txt,spades.log}
done

find -L "$1" -type d -regextype posix-extended \
 -regex '\/.*\/(K[1-9][0-9]{1,2}|misc|tmp)' -print | xargs rm -rv
