DEET 
[![DOI](https://zenodo.org/badge/19044/joshuaburkhart/DEET.svg)](https://zenodo.org/badge/latestdoi/19044/joshuaburkhart/DEET)
====

A program that finds and annotates genes using contigs and singletons yielded from mosquito microarray assays.

Input:

    1. A fasta file containing both contig & singleton sequences, specified with -f.

    2. A directory containing either Nimblegen Microarray or RNAseq data files, respectively specified with -m or -r.

Examples:

    1. ruby deet.rb -f ~/research/Bradshaw/gene_expression/microarray/fastas/Contigs_and_Singletons.simple.fasta -m ~/research/Bradshaw/gene_expression/microarray/ma_dat/photoperiod/

    2. ruby deet.rb -f ~/research/Bradshaw/gene_expression/rna_seq/fastas/contigs.fa -r ~/research/Bradshaw/gene_expression/rna_seq/rna_dat/

Algorithm Overview:

    1. parse provided fasta file(s), storing sequences longer than MIN_SEQ_LEN (set to 100 by default)
    2. report sequences shorter than or equal to MIN_SEQ_LEN in .seqs.csv
    2. hash expression data by sequence id, setting each value to that sequence's expression pattern (where 0=under expressed, 1=over expressed, and X=not significantly expressed)
    3. filter fasta sequences, assuring their existance in the hash
    4. query NCBI web interface for sequence
    5. record unmatched sequences and those with a low e value in .seqs.csv (E_LIM set to 1.0e-5 by default)
    6. group sequences by accession number, then by expression signature
    7. report groups of sequences in .result .seqs.csv files

