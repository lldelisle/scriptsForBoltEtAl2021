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

# Set the parameters:
directory <- "mcmc/"
directory.plot <- "./"
clusters.used <- c(0:9, 11)

# Here are the model which have been manually identified as not converging:
notconverging <- c()

# First we identify the files with the log evidence to identify the best model
logevid.files <- list.files(path = directory, pattern = "logevid.txt")
# We use the file names to get metadata
meta.data <- data.frame(matrix(unlist(strsplit(logevid.files, "_")), ncol = 5, byrow = T))
colnames(meta.data) <- c("model", "genotype", "cluster", "genes", "extension")
meta.data$logevid.file <- logevid.files
meta.data$extension <- NULL
# We set the order of the cluster and genotype
meta.data$cluster <- factor(meta.data$cluster, levels = paste0("cl", clusters.used))
meta.data$genotype <- factor(meta.data$genotype, levels = c("wt", "inv2"))

# We set the order of genes
meta.data <- cbind(meta.data, matrix(unlist(strsplit(meta.data$genes, "VS")), ncol = 2, byrow = T))
colnames(meta.data)[-1:0 + ncol(meta.data)] <- c("genex", "geney")
meta.data$geney <- factor(meta.data$geney, levels = c("Hoxd11", "Hoxa13", "Hoxd13"))
meta.data$genex <- factor(meta.data$genex, levels = c("Hoxa11", "Hoxd11", "Hoxa13"))

# We read the log evidence from the file
meta.data <- ddply(meta.data, .(logevid.file), mutate,
                   log.evid = as.numeric(readLines(file.path(directory, logevid.file), n = 1)))
# For the model which are not converging we put -Inf
meta.data$log.evid[gsub("_logevid.txt", "", meta.data$logevid.file) %in% notconverging] <- -Inf
# Identify the best model
meta.data <- ddply(meta.data, .(genotype, cluster, genes), mutate, is.best = log.evid == max(log.evid, na.rm = T))
meta.data$is.best[is.na(meta.data$is.best)] <- F
meta.data <- ddply(meta.data, .(genotype, cluster, genes), mutate, best = model[is.best])

# Check the best models
with(subset(meta.data, is.best), table(model))

# Create a new table with info condensed:
summary.log.evid.df <- cast(meta.data, formula = genotype + cluster + genes + best ~ model, value = "log.evid")

# At this step the convergence of best models is manually checked.
# If one model is not converging the notconverging variable is updated and 
# The lines above are re-run.

# The mean pdf is read of the best models
meta.data$file <- gsub("_logevid", "_pdf2d", meta.data$logevid.file)
pdfs <- do.call(rbind, lapply(meta.data$file[meta.data$is.best], function(fn){
  df <- read.delim(file.path(directory, fn))
  colnames(df)[1] <- "y"
  df2 <- melt(df, 
              id.vars = "y",
              variable_name = "x")
  df2$x <- as.numeric(gsub("X", "", df2$x))
  df2$file <- fn
  return(df2)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# Get the corr
meta.data$file.cor <- gsub("pdf2d", "corr", meta.data$file)
corr <- do.call(rbind, lapply(meta.data$file.cor[meta.data$is.best], function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file.cor <- fn
  return(df)
}))
# Add the meta data
corr <- merge(corr, meta.data)
# A label is formatted:
corr$label <- paste0(round(corr$mean, 2), "[-", round(corr$mean, 2) - round(corr$low, 2), "]",
                     "^{+", round(corr$high, 2) - round(corr$low, 2), "}")
# For the p-value, a superior value is given when significant
# a inferior value when not
corr$p.signif <- corr$pval + corr$error < 0.05
corr$p.n.signif <- corr$pval - corr$error > 0.05
with(corr, table(p.signif, p.n.signif))
corr$p.label <- sapply(corr$pval + corr$error, function(v){paste0("p<", format(v, digits = 2))})
corr$p.label[corr$p.n.signif] <- sapply(corr$pval[corr$p.n.signif] - corr$error[corr$p.n.signif],
                                        function(v){paste0("p>", format(v, digits = 2))})
corr$p.label[!corr$p.n.signif & !corr$p.signif] <- "p around 0.05"

# Figure 3f
# This is a multipanel figure:
my.df <- subset(pdfs, genex == "Hoxa13" & cluster == "cl3")
my.corr <- subset(corr, genex == "Hoxa13" & cluster == "cl3")
p1 <- ggplot(my.df, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  geom_text(data = my.corr,
            aes(label = label),
            x = 2.5, y = 2.5,
            size = 2.5,
            parse = T) +
  geom_text(data = my.corr,
            aes(label = p.label),
            x = 2.5, y = 1.8,
            size = 2.5) +
  facet_grid(. ~ genotype) +
  theme_minimal() +
  ylab("") +
  xlab(unique(my.df$genex)) +
  scale_fill_gradient(low="white", high="black", limits = c(0, 2.5)) + 
  theme(strip.text.x = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
        legend.position = "none",
        axis.title.x = element_text(face = 'italic'))

my.df <- subset(pdfs, genex == "Hoxd11" & cluster == "cl6")
my.corr <- subset(corr, genex == "Hoxd11" & cluster == "cl6")
p2 <- ggplot(my.df, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  geom_text(data = my.corr,
            aes(label = label),
            x = 2.5, y = 2.5,
            size = 2.5,
            parse = T) +
  geom_text(data = my.corr,
            aes(label = p.label),
            x = 2.5, y = 1.8,
            size = 2.5) +
  facet_grid(. ~ genotype) +
  theme_minimal() +
  ylab("") +
  xlab(unique(my.df$genex)) +
  scale_fill_gradient(low="white", high="black", limits = c(0, 2.5)) + 
  theme(strip.text.x = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
        legend.position = "none",
        axis.title.x = element_text(face = 'italic'))

my.df <- subset(pdfs, genex == "Hoxa11" & cluster == "cl3")
my.corr <- subset(corr, genex == "Hoxa11" & cluster == "cl3")
p3 <- ggplot(my.df, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  geom_text(data = my.corr,
            aes(label = label),
            x = 2.5, y = 2.5,
            size = 2.5,
            parse = T) +
  geom_text(data = my.corr,
            aes(label = p.label),
            x = 2.5, y = 1.8,
            size = 2.5) +
  facet_grid(. ~ genotype) +
  theme_minimal() +
  ylab("") +
  xlab(unique(my.df$genex)) +
  scale_fill_gradient(low="white", high="black", limits = c(0, 2.5)) + 
  theme(strip.text.x = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
        legend.position = "none",
        axis.title.x = element_text(face = 'italic'))

my.df <- subset(pdfs, genex == "Hoxa11" & cluster == "cl6")
my.corr <- subset(corr, genex == "Hoxa11" & cluster == "cl6")
p4 <- ggplot(my.df, aes(x, y)) +
  geom_tile(aes(fill = log(1 + value))) +
  geom_text(data = my.corr,
            aes(label = label),
            x = 2.5, y = 2.5,
            size = 2.5,
            parse = T) +
  geom_text(data = my.corr,
            aes(label = p.label),
            x = 2.5, y = 1.8,
            size = 2.5) +
  facet_grid(. ~ genotype) +
  theme_minimal() +
  ylab("") +
  xlab(unique(my.df$genex)) +
  scale_fill_gradient(low="white", high="black", limits = c(0, 2.5)) + 
  theme(strip.text.x = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ), 
        legend.position = "none",
        axis.title.x = element_text(face = 'italic'))
fig <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggexport(fig, filename = file.path(directory.plot, "fig3f.pdf"), width = 6.5, height = 3.7)
