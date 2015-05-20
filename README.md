#DEET
A program that finds and annotates genes using contigs and singletons yielded from mosquito microarray assays.

##Input
    1. A fasta file containing both contig & singleton sequences, specified with -f.

    2. A directory containing either Nimblegen Microarray or RNAseq data files, respectively specified with -m or -r.

##Examples
    $ ruby deet.rb -f ~/research/Bradshaw/gene_expression/microarray/fastas/Contigs_and_Singletons.simple.fasta -m ~/research/Bradshaw/gene_expression/microarray/ma_dat/photoperiod/

    $ ruby deet.rb -f ~/research/Bradshaw/gene_expression/rna_seq/fastas/contigs.fa -r ~/research/Bradshaw/gene_expression/rna_seq/rna_dat/

##Algorithm Overview
    1. parse provided fasta file(s), storing sequences longer than MIN_SEQ_LEN (set to 100 by default)
    2. report sequences shorter than or equal to MIN_SEQ_LEN in .seqs.csv
    2. hash expression data by sequence id, setting each value to that sequence's expression pattern (where 0=under expressed, 1=over expressed, and X=not significantly expressed)
    3. filter fasta sequences, assuring their existance in the hash
    4. query NCBI web interface for sequence
    5. record unmatched sequences and those with a low e value in .seqs.csv (E_LIM set to 1.0e-5 by default)
    6. group sequences by accession number, then by expression signature
    7. report groups of sequences in .result .seqs.csv files

##DEET Local DB Prep
###Reference
    http://www.ncbi.nlm.nih.gov/books/NBK1763/

###Download and Combine Files
    $ wget -r -nc -A "*.rna.fna.gz" -w 30 ftp://ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate/
    $ cd /N/u/joshburk/Mason/refseq_complete/ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate
    $ gunzip ./*
    $ cat ./* > invertebrate.rna.fna
    
###Make Blast DB
    https://www.biostars.org/p/97829/

    $ /N/u/joshburk/Mason/blast/ncbi-blast-2.2.29+-src/c++/ReleaseMT/bin/makeblastdb -in invertebrate.rna.fna -dbtype nucl -parse_seqids -out /N/u/joshburk/Mason/refseq_complete/invertebrate_rna.fna

###Make Header File
    $ cat invertebrate.rna.fna | grep -e '^>' > /N/u/joshburk/Mason/refseq_complete/invertebrate_rna.fna.headers_only

###Run DEET Local DB
    $ module load ruby
    $ ruby deet_local_db.rb \
    -f /N/u/joshburk/Mason/gene_expression/microarray/fastas/Contigs_and_Singletons.simple.fasta \
    -m /N/u/joshburk/Mason/gene_expression/new_photoperiod_ma-2015-05-05/new_photoperiod_ma/

###*run-deet_local.pbs*
    --
    #PBS -N deet_local
    #PBS -l walltime=128:00:00
    #PBS -q batch
    #PBS -l vmem=500gb
    #PBS -M burkhart.joshua@gmail.com
    #PBS -m abe
    #PBS -l nodes=1:ppn=32

    module load ruby

    cd ~/DEET && \
    ruby deet_local_db.rb \
    -f /N/u/joshburk/Mason/gene_expression/microarray/fastas/Contigs_and_Singletons.simple.fasta \
    -m /N/u/joshburk/Mason/gene_expression/new_photoperiod_ma-2015-05-05/new_photoperiod_ma/
    --
