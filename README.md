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

    1. parse simple fasta file, creating sequence objects for sequences 20bp or longer
    2. perform an ncbi blast with each sequence, storing top two alignments in a blast result object
    3. keep blast result objects for unique top alignment accession numbers (genes)
    4. keep blast result objects with top alignment e values < x (where -20 < x < -5)
    5. get sequence for top alignment accession number from NCBI
    6. submit sequence to AmiGO, linking through gene product info to uniprot to complete go annotation
    7. collect GO numbers and put into GODataRow objects
    8. create MADataRow objects from MA data
    9. create DoubleDataRow objects from GODataRow and MADataRow objects


