# Misc Genomics Scripts

- **calc.ATCG.content.bash**: a unix way to quickly get A,T,C,G,N,- content from a FastA file; handles linewraps and multiple records

- **filter.contigs.py**: cleans up a _de novo_ assembly from SPAdes, Velvet, or IDBA (requires biopython). IDBA includes space-delimited data in their contig headers, and because SeqIO parses on whitespace, these will need to be removed or replaced (e.g., `sed -i 's/ /|/g' assembly.fna`). SPAdes and Velvet lack whitespace in their contig deflines, so those output files can be directly fed into this filtering script.

    **Example batch usage**:
    If there are many SPAdes assembly output directories beginning with 3009 that need filtering, first cd into the parent dir containing all of the 3009\* dirs and execute: `for F in 3009*/contigs.fasta; do B=$(dirname $F | awk -F \/ '{print $1}'); filter.contigs.py -i $F -g -m -c 5 -l 250 -o "$B".fna; done` This took me 1 min for 340 assemblies.

- **find.dupes.bash**: given a FastA file, identifies repetitive regions, and outputs a BED file; flexible opts for defining repetitive sites; depends on BEDTools and MUMmer (nucmer)


***
# Fetching data from NCBI

- **biosample2FastQ.py**: download FastQ read files from NCBI for a given BioSample (or SRR) accession

- **genbankacc2gbk.py** A GenBank file is fetched from NCBI given an accession number. When more than one accession is provided (e.g., for taxa with >1 chromosome or harboring plasmids) all records are merged into a single output file. The `--min-length` option ensures unusually small sequence sizes don't make their way into downstream analyses. Biopython is used to avoid a dependency on efetch.

***
# GenBank file manipulations

- **extract.nucl.from.GBK.py**: Specify a search term such as a gene name or locus_tag and extract its nucleotide sequence. An ERROR message is printed if the query returns more than one or no hit.
- **gbk2proteome.molec.weights.py**: prints molecular weights of all proteins from coding sequences to stdout and provides the corresponding locus_tag and product for each as well. Handles unknown residues in proteins ('X') by estimating each as 128.16 Daltons and appends an '~' in front of the calculated molecular weight to indicate approximation.
- **locus_tag2faa.py**: given a locus tag, an amino acid FastA is printed to stdout

***
# Cleaning up disk space

- BLAST searching requires index files that can be easily and quickly re-generated, so remove all leftover binary files within $HOME: `find $HOME -type f -regextype posix-extended -regex '.*\.(nhr|nih|nin|nog|nsd|nsi|nsq|psi|psq)' -print | xargs rm -v`
- SPAdes keeps a lot of intermediate files, so delete these but keep essential log and FastA files to repeat the assembly if necessary: `prune.SPAdes.assembly.dirs.bash $HOME`
