options(stringsAsFactors = F)

# Install dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("ggrepel")
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# One of the dependency of ggpubr is only available for R=4.
if ( ! require("ggpubr") & ! rversionAbove(4)){
  packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  safelyLoadAPackageInCRANorBioconductor("ggpubr")
}
safelyLoadAPackageInCRANorBioconductor("ggrastr")

wd <- commandArgs(TRUE)[1]
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# The table is:
# input\tgenex\tgeney\txmax\tymax\tgroup
table.fn <- commandArgs(TRUE)[3]
# output prefix
output.prefix <- commandArgs(TRUE)[4]


setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "genex", "geney", "xmax", "ymax", "group"))
my.table$i <- rownames(my.table)

my.table$input <- basename(my.table$full.path.input)

my.table$genes <- paste0(my.table$genex, "VS", my.table$geney)

# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pdf2d_flat.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pdf2d_flat.txt", "", pdf.files), function(v){
  data <- strsplit(v, "_")[[1]]
  return(c(data[1:2], paste(data[-(1:2)], collapse = "_")))
})))
colnames(meta.data) <- c("model", "genes", "info")

meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)

meta.data$file <- pdf.files

meta.data$i <- sapply(strsplit(meta.data$info, "_"), tail, 1)

# Check the gene were correctly identified
temp <- merge(unique(my.table[, c("i", "genes")]), unique(meta.data[, c("i", "genes")]), by = "i")
if(!all(temp$genes.x == temp$genes.y)){
  # The gene had "_" in its name:
  meta.data$genes <- my.table$genes[match(meta.data$i, my.table$i)]
  meta.data$info <- apply(meta.data[, c("file", "genes")], 1, function(v){
    return(paste(strsplit(gsub("_pdf2d_flat.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$genes, "__", meta.data$info)
}
# Process genes
meta.data <- cbind(meta.data, matrix(unlist(strsplit(meta.data$genes, "VS")), ncol = 2, byrow = T))
colnames(meta.data)[-1:0 + ncol(meta.data)] <- c("genex", "geney")
# Find the group name
meta.data$group <- apply(meta.data[, c("id", "i")], 1, function(v){
  group.name <- my.table$group[my.table$i == v[2]]
  if (! is.na(group.name) & group.name != ""){
    return(paste0(group.name, gsub(paste0("_", v[2], "$"), "", strsplit(v[1], paste0("_", group.name))[[1]][2])))
  } else {
    return("all")
  }
})
meta.data$input <-my.table$input[match(meta.data$i, my.table$i)]
# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$value <- df$mean
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# Get the corr
meta.data$file.cor <- gsub("pdf2d_flat", "corr", meta.data$file)
corr <- do.call(rbind, lapply(meta.data$file.cor, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file.cor <- fn
  return(df)
}))
# Add the meta data
corr <- merge(corr, meta.data)
# A label is formatted:
corr$label <- paste0(round(corr$mean, 2), "[-", round(corr$mean, 2) - round(corr$low, 2), "]",
                     "^{+", round(corr$high, 2) - round(corr$mean, 2), "}")
# For the p-value, a superior value is given (still only mean + sd)
corr$p.label <- sapply(corr$pval + corr$error, function(v){paste0("p<", format(v, digits = 2))})

# I put wt first
my.groups <- c(grep("wt", unique(meta.data$group), ignore.case = T, value = T),
               grep("wt", unique(meta.data$group), ignore.case = T, value = T, invert = T))
my.labeller <- gsub("genotype", "", my.groups)
names(my.labeller) <- my.groups
pdfs$group <- factor(pdfs$group, levels = my.groups)
corr$group <- factor(corr$group, levels = my.groups)
ggplot.list <- list()
current <- 1
for (my.i in c(2, 4, 1, 3)){
  my.df <- subset(pdfs, i == my.i)
  my.corr <- subset(corr, i == my.i)
  ggplot.list[[current]] <- ggplot(my.df, aes(x, y)) +
    geom_tile_rast(aes(fill = log(1 + value))) +
    geom_text(data = my.corr,
              aes(label = label),
              x = 2.5, y = 2.5,
              size = 2.5,
              parse = T) +
    geom_text(data = my.corr,
              aes(label = p.label),
              x = 2.5, y = 1.8,
              size = 2.5) +
    facet_grid(input ~  group, labeller = labeller(group = my.labeller)) +
    theme_minimal() +
    ylab("") +
    xlab(unique(my.df$genex)) +
    scale_fill_gradient(low="white", high="black", limits = c(0, 2.5)) + 
    theme(# strip.text.x = element_blank(),
          plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
          legend.position = "none",
          axis.title.x = element_text(face = 'italic'))
  current <- current + 1
  
}
ggarrange(plotlist = ggplot.list)
ggplot.list.clean <- lapply(ggplot.list, function(gg){
  gg +
    theme(strip.text.x = element_blank(),
          strip.text.y = element_blank())
})
g <- ggarrange(plotlist = ggplot.list.clean)

ggsave(paste0(output.prefix, ".pdf"), g,  width = 6.5, height = 3.7)
