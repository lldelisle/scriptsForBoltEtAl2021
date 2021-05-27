options(stringsAsFactors = F)

# Install dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("ggrepel")

delog <- function(x) (exp(x) - 1) * 1e-4


wd <- commandArgs(TRUE)[1]
# The inputs are in the directory
directory <- commandArgs(TRUE)[2]
# The table is:
# input\tgene\txmax\tgroup
table.fn <- commandArgs(TRUE)[3]
# output prefix
output.prefix <- commandArgs(TRUE)[4]


setwd(wd)

# Get the table used as input of the mcmc
my.table <- read.delim(table.fn, h=F,
                       col.names = c("full.path.input", "gene", "xmax", "group"))
my.table$i <- rownames(my.table)
my.table$input <- basename(my.table$full.path.input)
# Find pdf files
pdf.files <- list.files(path = directory, pattern = "1-.*pdf.txt")

# Analyze pdf file name to get meta data
meta.data <- data.frame(t(sapply(gsub("_pdf.txt", "", pdf.files), function(v){
  data <- strsplit(v, "_")[[1]]
  return(c(data[1:2], paste(data[-(1:2)], collapse = "_")))
})))
colnames(meta.data) <- c("model", "gene", "info")

meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)

meta.data$file <- pdf.files

meta.data$i <- sapply(strsplit(meta.data$info, "_"), tail, 1)

# Check the gene were correctly identified
temp <- merge(unique(my.table[, c("i", "gene")]), unique(meta.data[, c("i", "gene")]), by = "i")
if(!all(temp$gene.x == temp$gene.y)){
  # The gene had "_" in its name:
  meta.data$gene <- my.table$gene[match(meta.data$i, my.table$i)]
  meta.data$info <- apply(meta.data[, c("logevid.file", "gene")], 1, function(v){
    return(paste(strsplit(gsub("_pdf.txt", "", v[1]), paste0(v[2], "_"))[[1]][-1], collapse = "_"))
  })
  meta.data$id <- paste0(meta.data$gene, "__", meta.data$info)
}
# I put the order as in my.table
my.table$gene <- factor(my.table$gene, levels = unique(my.table$gene))
meta.data$gene <- factor(meta.data$gene, levels = levels(my.table$gene))
# Find the group name
meta.data$group <- apply(meta.data[, c("id", "i")], 1, function(v){
  group.name <- my.table$group[my.table$i == v[2]]
  if (! is.na(group.name) & group.name != ""){
    return(paste0(group.name, gsub(paste0("_", v[2], "$"), "", strsplit(v[1], paste0("_", group.name))[[1]][2])))
  } else {
    return("all")
  }
})
# Add the input
meta.data$input <- my.table$input[match(meta.data$i, my.table$i)]

# Read the pdfs
pdfs <- do.call(rbind, lapply(meta.data$file, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# We evaluate FC and means:
fc.results <- NULL
means.results <- NULL

# I will use a 68% confidence interval
lvl <- 0.6827
qm <- (1 - lvl) / 2
qp <- 1 - qm

my.groups <- c(grep("wt", unique(meta.data$group), ignore.case = T, value = T),
               grep("wt", unique(meta.data$group), ignore.case = T, value = T, invert = T))
set.seed(1)
for (i in 1:nrow(my.table)){
  print(i)
  my.groups.i <- intersect(my.groups, meta.data$group[meta.data$i == i])
  if (length(my.groups.i) == 2){
    # I get the data
    means <- lapply(my.groups.i, function(g){
      as.numeric(readLines(file.path(directory, gsub("pdf.txt", "means.txt.gz",
                                                     meta.data$file[meta.data$i == i & meta.data$group == g]))))
    })
    # I get the size
    min.size <- min(sapply(means, length))
    # Subset
    means <- lapply(means, function(v){
      v[round(seq(1, length(v), length.out = min.size))]
    })
    means <- as.data.frame(means)
    colnames(means) <- my.groups.i
    means[, 2] <- sample(means[, 2], nrow(means), replace = F)
    means$logFC <- means[, 2] - means[, 1]
    means$FC <- delog(means[, 2]) / delog(means[, 1])
    means.results <- rbind(means.results,
                           c(i, colnames(means)[1],
                             mean(means[, 1]), # Average of log expression first group
                             quantile(means[, 1], probs = c(qm, 0.5, qp)), # Quantile on log expression first group
                             mean(delog(means[, 1])), # Average of delog expression first group ! this might be bad if extreme values exists
                             quantile(delog(means[, 1]), probs = c(qm, 0.5, qp)) # Quantile on log expression first group
                           ),
                           c(i, colnames(means)[2],
                             mean(means[, 2]), # Average of log expression second group
                             quantile(means[, 2], probs = c(qm, 0.5, qp)), # Quantile on log expression second group
                             mean(delog(means[, 2])), # Average of delog expression second group ! this might be bad if extreme values exists
                             quantile(delog(means[, 2]), probs = c(qm, 0.5, qp)) # Quantile on log expression second group
                           ))
    fc.results <- rbind(fc.results,
                        c(i, my.groups.i,
                          mean(means$logFC), # Average of log expression FC
                          quantile(means$logFC, probs = c(qm, 0.5, qp)), # Quantile on log FC
                          mean(means$FC), # Average of delog expression FC ! this might be bad if extreme values
                          quantile(means$FC, probs = c(qm, 0.5, qp)), # Quantile on delog FC
                          sum(means$FC > 1) / nrow(means), # Prop increased
                          sum(means$FC < 1) / nrow(means) # Prop decreased
                        ))
  }
}
fc.results <- data.frame(fc.results)
colnames(fc.results) <- c("i", "group1", "group2",
                          "meanlogFC", "lowlogFC", "medlogFC", "highlogFC",
                          "meanFC", "lowFC", "medFC", "highFC",
                          "propIncr", "propDecr")
fc.results <- merge(fc.results, my.table)
means.results <- as.data.frame(means.results)
colnames(means.results) <- c("i", "group",
                             "meanlog", "lowlog", "medlog", "highlog",
                             "mean", "low", "med", "high")
means.results <- merge(means.results, meta.data[, c("i", "group", "input", "gene")])

fc.results$value.label <- with(fc.results, paste0(format(round(as.numeric(medFC), 2), nsmall = 2),
                                                  "[", format(round(as.numeric(lowFC), 2), nsmall = 2), "]",
                                                  "^{", format(round(as.numeric(highFC), 2), nsmall = 2), "}"))
input.labeller <- gsub(".txt", "", gsub("_", " ", unique(my.table$input)))
names(input.labeller) <- unique(my.table$input)
pdfs$group <- factor(pdfs$group, levels = my.groups)
g <- ggplot(pdfs, aes(x = x, color = group, fill = group)) +
  geom_path(aes(y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.3, color = NA) +
  facet_grid(input ~ gene, labeller = labeller(input = input.labeller)) +
  geom_text(data = fc.results,
            inherit.aes = F,
            aes(label = value.label),
            x = 2.5, y = 2.5, parse = T) +
  theme_classic() +
  theme(legend.position = "top") +
  xlab("log(1+10'000*expression)") + 
  coord_cartesian(ylim=c(0, 3)) +
  ylab("Density") +
  scale_color_discrete("genotype",
                     labels = gsub("genotype", "", my.groups)) +
  scale_fill_discrete("genotype",
                       labels = gsub("genotype", "", my.groups))
  
ggsave(paste0(output.prefix, ".pdf"), g, width = 8, height = 4.5)
