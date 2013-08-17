deet
====

A program that finds and annotates genes using contigs and singletons yielded from mosquito microarray assays.

Input format:

    A fasta file containing both contig & singleton sequences with a ">" and a sequence identifier on one line and the corresponding sequence on the next line. 

    Example:

    ...

    >CONTIG06215
    ACAATTCTAGTAAAGCATTTGTCATCTTCTTTATATTTTTGGCAGTAaaCAGTCTGCAGTTCACCGCTCAACTCTTCGGCGCAACTTCTGTTCGTATCATTAGTTTTACATTGGTAACACCGGAGACGATCCATTGGGAAGACATCCCTATTGCAGTTGTCTGTTCCACAAACGCTGCAGGTTTTACCTGTACAGGCTTGCGTCTCGTTCAGCGTTAGATTTCCTTCACAGCCTCGAATGACATGGCCATCATCAACTACTTTCGAGTAACATTTCCCTGCAAACGCTTTAGCGCACGGTCCACTGCTGGTTGTGCCCAGAAGACAcTTTtCGTCATCGTCGTCGTCAACTTCCGAACTACATTGCAAACATTGCCGGACAGCCGGAGCCTCTCGATTGCAACCGTCAGTAGAGCACGTtCTGcATTCCCTACTATCCagACAAGCATCCACATTCAAACCAAGATCATTCTCACACCCTCTAGTGAGAgttGCACCCTCCATCCGTTCGTAGCAtcGGTTGTTACTcTTCCAATACTGACAGTAACCGGTGGCGTTGGTTTGCTGCTCagCACACGTGGCATCGGAAGCCTGGTTGCACTTGTAACACTGTAACCAATGGTCGGTATTACAGCCAGCAGAGTCACAGACAGTACAAGTTGCGTCCTCGGTTCGATTGCAAAGTGTTTGCTCTACAGCCGGAAGTGTTGACAAGCAACCGCGCTCAAGAAGGCCTCCTGTAACTCGTGTATAGCAGTTGTCCTGCTGTGTGATACATTTCTTCGGTGCGGCTGTACCATCCACACAGGTAGCGTCGGTTGAATTACACTGCAGACACCGGGAAGAGTTCTTCATCTCGGtgACGGACAGcTtGTTGCACAaTTCACCCGAGCAAATCGTACACACATCATTGTCACACGCGGTTTGCGGATCGGCGAAATCCGACAAACAGCCTCTGGTAACATTATCGCCTTCTATTTTAGTGTAACATCGATCATCGCTGCTGAAAATCGAACAAAGCTTGCCTCCCcTTTCCGGTTGCTCGGTGgTACAGTTTGCAGTGTTGGTCGAATCGCAGCTATGACACTTCAACCAGGGTGCGTTGTTACATCCGGCCCcGGTAcAGGTTACGCAGCTTGAATCAGTGCTATtGGTACAAGCACTAGCGATCGGATCcACCAGGCAGCCTCGTACGACGGCGTTGTCACTGGTTAGCATAGTGTAGCATTGATTAGCC
    >F5BTJ3O01C7KYF
    AATCACGACCAGCCAAGTCCAGACGGAGGATGGGCATGTGGCAAAGGGCACACAGGGG
    >F5BTJ3O01BEOZD
    CAACATCCGCAAGCTGTACAATCTAGGCAAGGATGACGATGTGCGTCAGTTCGTCGTTAA

    ...

    Procedure using https://github.com/joshuaburkhart/bio programs:

    $fa2oneline.pl Contigs.fna > Contigs.oneline.fna
    $fa2oneline.pl Singletons.fna > Singletons.oneline.fna
    $fasta_cleaner.sh Contigs.oneline.fna 6
    $fasta_cleaner.sh Singletons.oneline.fna 7
    $cat Contigs.oneline.fna.indexed | awk -F' ' '{print ">"$1"\n"$2}' > Contigs.simple.fasta
    $cat Singletons.oneline.fna.indexed | awk -F' ' '{print ">"$1"\n"$2}' > Singletons.simple.fasta
    $cp Contigs.simple.fasta Contigs_and_Singletons.simple.fasta
    $cat Singletons.simple.fasta >> Contigs_and_Singletons.simple.fasta

Algorithm Overview:

    1. work through Justin Choi's pipeline (doi: 10.1186/1471-2164-11-703)
    2. parse Justin Choi's output along with FASTA files, representing clusters by their longest sequence
    3. tblastx each sequence against Arthropods, Anopheles, and Drosophila
    4. record the two alignments for each ncbi result with the lowest E values
    5. group ncbi results by accession number of alignment with lowest E value
    6. parse microarray data, dividing accession number groups by expression signature
    7. use accession number groups from 5 and expression signature groups from 6 to produce a list of putative genes and paralogs
    8. represent expression signature groups by their longest sequence
    9. produce 2D expression plots
    10. submit sequences to U. Indiana's annotation web tool
    11. submit resulting annotations to a pie chart tool


