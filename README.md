# Misc Genomics Scripts
- **calc_ATCG_content.bash**: a unix way to quickly get A,T,C,G,N,- content from a FastA file; handles linewraps and multiple records

- **filter_contigs.py**: cleans up a _de novo_ assembly from SPAdes, Velvet, or IDBA (requires biopython). IDBA includes space-delimited data in their contig headers, and because SeqIO parses on whitespace, these will need to be removed or replaced (e.g., `sed -i 's/ /|/g' assembly.fna`). SPAdes and Velvet lack whitespace in their contig deflines, so those output files can be directly fed into this filtering script.

    **Example batch usage**:
    If there are many SPAdes assembly output directories beginning with 3009 that need filtering, first cd into the parent dir containing all of the 3009\* dirs and execute: `for F in 3009*/contigs.fasta; do B=$(dirname $F | awk -F \/ '{print $1}'); filter_contigs.py -i $F -g -m -c 5 -l 500 -o "$B".fna -p ~/SPAdes_Assems/filtered_contigs; done` This took me 1 min for 340 assemblies.

- **find_dupes.bash**: given a FastA file, identifies repetitive regions, and outputs a BED file; flexible opts for defining repetitive sites; depends on BEDTools and MUMmer (nucmer)

***

# Cleaning up disk space
- BLAST searching requires index files that can be easily and quickly re-generated, so remove all leftover binary files within the Desktop dir: `find ~/Desktop -type f -regextype posix-extended -regex '.*\.(nih|nin|nsq|psi|psq)' -print | xargs rm -v`
- SPAdes keeps a lot of intermediate files, so delete these but keep essential log and FastA files to repeat the assembly if necessary: `bash ~/genomics_scripts/prune_SPAdes_assembly_dirs.bash ~/Desktop`
