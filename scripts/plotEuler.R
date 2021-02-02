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
safelyLoadAPackageInCRANorBioconductor("eulerr")
if(length(commandArgs(TRUE))>0){
  if (length(commandArgs(TRUE)) < 2){
    cat("Usage: Rscript plotEuler.R fileWithNb.txt outputPDF")
    stop()
  }
  f <- commandArgs(TRUE)[1]
  output <- commandArgs(TRUE)[2]
} else{
  #Ask for the config file
  print("Provide the file with line counts.")
  f <- file.choose()
  output <- file.path(dirname(f), paste0(file_path_sans_ext(f), ".pdf"))
}
df <- read.table(f)
df <- subset(df, V2 != "total")
df$expe <- basename(df$V2)
single.names <- setdiff(grep("vs", df$expe, invert = T, value = T), df$expe[nrow(df)])
df$expe <- gsub("vs", "&", df$expe)
df$expe[nrow(df)] <- paste(single.names, collapse = "&")
quantities <- df$V1
names(quantities) <- df$expe
fit.euler <- euler(quantities, input = "union")
meta.data <- as.data.frame(matrix(unlist(strsplit(single.names, "_")), ncol = 5, byrow = T))
meta.data$V5 <- gsub(".bed", "", meta.data$V5)
pretty.names <- paste(meta.data$V5,
                      meta.data$V3,
                      format(df$V1[match(single.names, df$expe)], big.mark = ","),
                      sep = "\n")
pdf(output)
plot(fit.euler, fill = c("blue", "grey", "purple"), labels = pretty.names, quantities = T,)
dev.off()
