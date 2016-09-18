
# DEET Post Processing R Script
# Author: Joshua Burkhart
# Date: 9/17/2016

# Libraries
library(openxlsx)
library(dplyr)
library(magrittr)

# Globals
LMMA_DIR <- "/Users/joshuaburkhart/Research/DEET/photoperiod/new_ma/"
SEQS_CSV_DIR <- "/Users/joshuaburkhart/Research/DEET/photoperiod/"
SEQS_CSV_FILE <- "1460605014.seqs.csv"
SEQS_CSV_PATH <- paste(SEQS_CSV_DIR,SEQS_CSV_FILE,sep="")
OUT_DIR <- SEQS_CSV_DIR
OUT_FILE <- paste(SEQS_CSV_FILE,".post_processed.xlsx",sep="")
OUT_PATH <- paste(OUT_DIR,OUT_FILE,sep="")

# Read seqs.csv into dataframe (sep="~")
seqs.csv <- read.delim(SEQS_CSV_PATH,header=TRUE,sep="~")

setwd(LMMA_DIR)
temp = list.files(pattern="*.txt") # bc list.files doesn't actually work

# Read all limma files in dir into dataframes (sep="\t")
# Append limma M, A, pvalue, qvalue columns to seqs.csv, matching gene to ID
# Label columns like "limma.xxx.ABC-DEF.gene.de.txt M" 
for (i in 1:length(temp)) {
  assign(temp[i], read.delim(temp[i],header=TRUE,sep="\t"))
  
  # resolve string name
  cur.limma <- get(temp[i])
  
  # assure matching case
  cur.limma$gene <- cur.limma$gene %>% toupper()
  seqs.csv$ID <- seqs.csv$ID %>% toupper()
  
  # rename cols of interest
  names(cur.limma)[c(2,3,5,6)] <- paste(temp[i],names(cur.limma)[c(2,3,5,6)]) # only M, A, pvalue, qvalue column idx's
  
  # select id column and cols of interest
  cur.limma <- cur.limma[c(1,2,3,5,6)]
  
  # join
  seqs.csv <- seqs.csv %>% dplyr::left_join(cur.limma,by = c("ID"="gene"))
}

# Write output to xlsx
seqs.csv %>% openxlsx::write.xlsx(file=OUT_PATH)
