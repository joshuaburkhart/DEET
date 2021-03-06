

library(dplyr)
library(rBLAST)
library(seqinr)
library(assertthat)

TRANSCRIPTOME_DIR <-
  "/media/burkhart/Media/Research/Bradshaw/DEET/transcriptomes/"
QRY_TRANSCRIPTOME <-
  paste(TRANSCRIPTOME_DIR, "454wyeomyia.CONTIGnATCG.fa", sep = "")
RF1_TRANSCRIPTOME <-
  paste(TRANSCRIPTOME_DIR, "Aedes_aegypti.AaegL3.cds.all.fa", sep = "")
RF2_TRANSCRIPTOME <-
  paste(TRANSCRIPTOME_DIR,
        "Anopheles_gambiae.AgamP4.cds.all.fa",
        sep = "")
RF3_TRANSCRIPTOME <-
  paste(TRANSCRIPTOME_DIR,
        "Culex_quinquefasciatus.CpipJ2.cds.all.fa",
        sep = "")
# generated with:
# https://github.com/joshuaburkhart/bio
#$ ma_inverter.rb biting/ma/limma.WIOB-KC.gene.de.txt biting/ma/limma.KC-WIOB.gene.de.txt
#$ ma2plottable_csv.rb biting/ma/limma.WI-WIOB.gene.de.txt biting/ma/limma.KC-WIOB.gene.de.txt biting/analysis/1450736707.seqs.csv
QUAD_DATA <-
  "/home/burkhart/Research/Bradshaw/DEET/WI-WIOB-KC-WIOB.csv"

# first spreadsheet from workbook 4
QRY_SEQ_IDS <-
  "/home/burkhart/Research/Bradshaw/DEET/QuadPPlotGenes from John 29 Nov 2016 IV odorant.csv"
SINGLETON_REGEX <- "^F5BTJ3O0"
CONTIG_REGEX <- "^CONTIG"
OUTPUT_DIR <- "/home/burkhart/Software/DEET/output/"
RESULTS_FILEPATH <- paste(OUTPUT_DIR,"quadPlotHomologueAlignments.csv",sep="")

# read quad data
quad_df <- read.csv(QUAD_DATA,
                    header = TRUE,
                    stringsAsFactors = FALSE)

# read query seq ids using regex
qry_seq_ids_df <- read.csv(QRY_SEQ_IDS,
                           header = TRUE,
                           stringsAsFactors = FALSE) %>%
  dplyr::filter(
    grepl(SINGLETON_REGEX, QuadPlotGenes.from.John.29.Nov.2016.II) |
      grepl(CONTIG_REGEX, QuadPlotGenes.from.John.29.Nov.2016.II)
  )

# read transcriptomes
qry_transcriptome <- seqinr::read.fasta(QRY_TRANSCRIPTOME)
rf1_transcriptome <- seqinr::read.fasta(RF1_TRANSCRIPTOME)
rf2_transcriptome <- seqinr::read.fasta(RF2_TRANSCRIPTOME)
rf3_transcriptome <- seqinr::read.fasta(RF3_TRANSCRIPTOME)

# filter query transcriptome by seq ids
qry_transcriptome_filt <-
  qry_transcriptome[names(qry_transcriptome) %in%
                      qry_seq_ids_df$QuadPlotGenes.from.John.29.Nov.2016.II]

# create and load reference blast databases
rBLAST::makeblastdb(RF1_TRANSCRIPTOME)
rf1_db <- rBLAST::blast(RF1_TRANSCRIPTOME, type = "tblastx")

rBLAST::makeblastdb(RF2_TRANSCRIPTOME)
rf2_db <- rBLAST::blast(RF2_TRANSCRIPTOME, type = "tblastx")

rBLAST::makeblastdb(RF3_TRANSCRIPTOME)
rf3_db <- rBLAST::blast(RF3_TRANSCRIPTOME, type = "tblastx")

results_df <- data.frame()

itr <- 0
tot <- length(qry_transcriptome_filt)
for (qry_seq_id in names(qry_transcriptome_filt)) {
  
  # blast query seqs to references
  rf1_res <- predict(rf1_db,
                     Biostrings::DNAStringSet(seqinr::c2s(qry_transcriptome_filt[[qry_seq_id]])))
  
  assertthat::assert_that(nrow(rf1_res) >= 1)
  assertthat::assert_that(ncol(rf1_res) >= 11)
  
  # match
  rf1_match = rf1_res[1, 2]
  # e
  rf1_e = rf1_res[1, 11]
  
  rf2_res <- predict(rf2_db,
                     Biostrings::DNAStringSet(seqinr::c2s(qry_transcriptome_filt[[qry_seq_id]])))
  
  assertthat::assert_that(nrow(rf2_res) >= 1)
  assertthat::assert_that(ncol(rf2_res) >= 11)
  
  # match
  rf2_match = rf2_res[1, 2]
  # e
  rf2_e = rf2_res[1, 11]
  
  rf3_res <- predict(rf3_db,
                     Biostrings::DNAStringSet(seqinr::c2s(qry_transcriptome_filt[[qry_seq_id]])))
  
  assertthat::assert_that(nrow(rf3_res) >= 1)
  assertthat::assert_that(ncol(rf3_res) >= 11)
  
  # match
  rf3_match = rf3_res[1, 2]
  # e
  rf3_e = rf3_res[1, 11]
  
  assertthat::assert_that(qry_seq_id %in% quad_df$Sequence.ID)
  
  # x coord
  x_coord <-
    quad_df %>% dplyr::filter(Sequence.ID == qry_seq_id) %>% .[1, 2]
  # y coord
  y_coord <-
    quad_df %>% dplyr::filter(Sequence.ID == qry_seq_id) %>% .[1, 3]
  
  # create df with query seq, x coord, y coord, rf1 id, rf1 e, rf2 id, rf2 e, rf3 id, rf3 e
  results_df <- rbind(results_df,
                      data.frame(SequenceID = qry_seq_id,
                                 Aedes_aegypti_ID = rf1_match,
                                 Aedes_aegypti_e = rf1_e,
                                 Anopheles_gambiae_ID = rf2_match,
                                 Anopheles_gambiae_e = rf2_e,
                                 Culex_quinquefasciatus_ID = rf3_match,
                                 Culex_quinquefasciatus_e = rf3_e,
                                 quadplot_x_coordinate = x_coord,
                                 quadplot_y_coordinate = y_coord))
  
  if(itr %% 100 == 0){
    print(paste("processed ",itr," of ",tot," query sequences..."))
  }
  itr <- itr + 1
}

results_df <- results_df %>%
  dplyr::mutate(Quadrant = ifelse(
    quadplot_x_coordinate < 0,
    ifelse(quadplot_y_coordinate < 0,
           "LL gene",
           "UL gene"),
    ifelse(quadplot_y_coordinate < 0,
           "LR gene",
           "UR gene")
  ))

# save results to file
write.csv(results_df,RESULTS_FILEPATH)

# create example quadplot for spotchecking with quadrant counts
palette(c("firebrick2","black","black","deepskyblue"))
plot(x=results_df$quadplot_x_coordinate,
     y=results_df$quadplot_y_coordinate,
     pch="O",
     col=as.factor(results_df$Quadrant),
     xlab="Selection",
     ylab="End Points of Evolution",
     xlim=c(-6.9,6.9),
     ylim=c(-6.9,6.9))
abline(a=0,b=0)
abline(v=0)
#grid()
legend("topright",
       legend=sum(results_df$quadplot_x_coordinate > 0 &
                               results_df$quadplot_y_coordinate > 0),
       bty="n")
legend("topleft",
       legend=sum(results_df$quadplot_x_coordinate < 0 &
                              results_df$quadplot_y_coordinate > 0),
       bty="n")
legend("bottomleft",
       legend=sum(results_df$quadplot_x_coordinate < 0 &
                                 results_df$quadplot_y_coordinate < 0),
       bty="n")
legend("bottomright",
       legend=sum(results_df$quadplot_x_coordinate > 0 &
                                  results_df$quadplot_y_coordinate < 0),
       bty="n")

