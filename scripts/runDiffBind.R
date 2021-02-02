options(scipen=999)
options(stringsAsFactors=F)
rm(list=ls())
if(! "devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
library(tools)
#Install the libraries needed if they are not already present
safelyLoadAPackageInCRANorBioconductor("DiffBind")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
if(length(commandArgs(TRUE))>0){
  if (length(commandArgs(TRUE)) < 2){
    cat("Usage: Rscript runDiffBind.R sampleSheet.csv outputNoExtension")
    stop()
  }
  f <- commandArgs(TRUE)[1]
  output <- commandArgs(TRUE)[2]
} else{
  #Ask for the config file
  print("Provide the sampleSheet")
  f <- file.choose()
  output <- file.path(dirname(f), file_path_sans_ext(f))
}
# Load the samples:
samples_dba <- dba(sampleSheet=f, config=NULL)
# Count the reads:
samples_count <- dba.count(samples_dba, summits=250)
samples_contrast <- dba.contrast(samples_count, minMembers=2, categories=c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION))
samples_analyze <- dba.analyze(samples_contrast, method=DBA_ALL_METHODS, bFullLibrarySize=TRUE, bCorPlot=TRUE)
samples_DB_report <- dba.report(samples_analyze, contrast=1, method=DBA_DESEQ2, DataType=DBA_DATA_GRANGES)
# Export the up:
export.bed(samples_DB_report[samples_DB_report$Fold > 0],
           con=paste0(output, "_", samples_contrast$contrasts[[1]]$name1, ".bed"))
# Export the down:
export.bed(samples_DB_report[samples_DB_report$Fold < 0],
           con=paste0(output, "_", samples_contrast$contrasts[[1]]$name2, ".bed"))
samples_NotDB_report <- dba.report(samples_analyze, contrast=1, method=DBA_DESEQ2, bNotDB=TRUE,  DataType=DBA_DATA_GRANGES)
samples_NotDB_GR <- dba.peakset(samples_NotDB_report, bRetrieve=TRUE)
# Export the no change:
export.bed(samples_NotDB_GR,
           con=paste0(output, "_NDB.bed"))
