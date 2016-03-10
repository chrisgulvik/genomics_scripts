# Misc Genomics Scripts
- **filter_contigs.py**: cleans up a _de novo_ assembly from SPAdes, Velvet, or IDBA (requires biopython). IDBA includes space-delimited data in their contig headers, and because SeqIO parses on whitespace, these will need to be removed or replaced (e.g., `sed -i 's/ /|/g' assembly.fna`). SPAdes and Velvet lack whitespace in their contig deflines, so those output files can be directly fed into this filtering script.

    **Example batch usage**:
    If there are many SPAdes assembly output directories beginning with 3009 that need filtering, first cd into the parent dir containing all of the 3009\* dirs and execute: `for F in 3009*/contigs.fasta; do B=$(dirname $F | awk -F \/ '{print $1}'); filter_contigs.py -i $F -g -m -c 5 -l 500 -o "$B".fna -p ~/SPAdes_Assems/filtered_contigs; done` This took me 1 min for 340 assemblies.
