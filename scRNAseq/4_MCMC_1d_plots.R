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

# Set the parameters:
directory <- "mcmc/"
directory.plot <- "./"
clusters.used <- c(0:9, 11)

# Here are the model which have been manually identified as not converging:
notconverging <- c(# "2gauss_wt_cl1_Hoxd11", the xlim has been extended to 3.5 to make it converge.
  "3gauss_inv2_cl6_Hoxa11",
  "4gauss_inv2_cl6_Hoxa11")

delog <- function(x) (exp(x) - 1) * 1e-4

# First we identify the files with the pdf
pdf.files <- list.files(directory, pattern = "_pdf.txt")
meta.data <- data.frame(matrix(unlist(strsplit(pdf.files, "_")), ncol = 5, byrow = T))
# We use the file names to get metadata
colnames(meta.data) <- c("model", "genotype", "cluster", "gene", "extension")
meta.data$file <- pdf.files
meta.data$extension <- NULL
# We set the order of the cluster, gene and genotype
meta.data$cluster <- factor(meta.data$cluster, levels = paste0("cl", c(0:9, 11)))
meta.data$gene <- factor(meta.data$gene, levels = c("Hoxd13", "Hoxa11", "Hoxd11", "Hoxa13", "Shox2"))
meta.data$genotype <- factor(meta.data$genotype, levels = c("wt", "inv2"))
# We read the pdf values
pdfs <- do.call(rbind, lapply(pdf.files, function(fn){
  df <- read.delim(file.path(directory, fn))
  df$file <- fn
  return(df)
}))
# Add the meta data
pdfs <- merge(pdfs, meta.data)

# Add new variable used for plotting
pdfs$sign <- (as.integer(pdfs$genotype != "wt") * 2 - 1)
pdfs <- ddply(pdfs, .(file), mutate, delta.x = min(x[x != min(x)]) - min(x))

# Get log evidence value
log.evid.fn <- list.files(directory, pattern = "_logevid.txt")
log.evid <- as.numeric(sapply(log.evid.fn, function(fn){readLines(file.path(directory, fn))}))
log.evid.df <- data.frame(fn = log.evid.fn, log.evid = log.evid)

# Add the meta data
log.evid.df$file <- gsub("logevid", "pdf", log.evid.df$fn)
log.evid.df <- merge(log.evid.df, meta.data)

# For the model which are not converging we put -Inf
log.evid.df$log.evid[gsub("_logevid.txt", "", log.evid.df$fn) %in% notconverging] <- -Inf

# Look for the best model
log.evid.df <- ddply(log.evid.df, .(genotype, cluster, gene), mutate, is.best = log.evid == max(log.evid))
log.evid.df <- ddply(log.evid.df, .(genotype, cluster, gene), mutate, best = model[is.best])

# Check the best models
with(subset(log.evid.df, is.best), table(model))

# Create a new table with info condensed:
summary.log.evid.df <- cast(log.evid.df, formula = genotype + cluster + gene + best ~ model, value = "log.evid")

# I will evaluate the fold-change
# I do the mean.expression using the mean pdf only on the best model:
my.df <- subset(pdfs, file %in% log.evid.df$file[log.evid.df$is.best])
mean.expression <- ddply(my.df, .(model, genotype, cluster, gene), summarize, meanlogexpr = sum(mean * x * delta.x))
mean.expression <- ddply(mean.expression, .(model, genotype, cluster, gene), mutate, meanexpr = delog(meanlogexpr))

# I select only the clusters where I have data for both genotypes
clusters <- unique(subset(mean.expression, select = c("genotype", "cluster")))
mean.expression.FC <- subset(mean.expression, cluster %in% clusters$cluster[duplicated(clusters$cluster)])

# I estimate the fold-change using the mean pdf
mean.expression.FC <- ddply(mean.expression.FC, .(cluster, gene), mutate, meanlogexprFC = meanlogexpr[genotype != 'wt'] - meanlogexpr[genotype == 'wt'])
mean.expression.FC <- ddply(mean.expression.FC, .(cluster, gene), mutate, meanexprFC = meanexpr[genotype != 'wt'] / meanexpr[genotype == 'wt'])
mean.expression.FC$file <- with(mean.expression.FC, paste(model, genotype, cluster, gene, "means.txt", sep = "_"))

# I will use a 68% confidence interval
lvl <- 0.6827
qm <- (1 - lvl) / 2
qp <- 1 - qm

# For the mean expression
# For the logFC (with pseudo count of 1)
# For the FC
mean.expression.FC$lowmeanexpr.shuffle <- NA
mean.expression.FC$meanexpr.shuffle <- NA
mean.expression.FC$medianexpr.shuffle <- NA
mean.expression.FC$highmeanexpr.shuffle <- NA
mean.expression.FC$lowlogexprFC.shuffle <- NA
mean.expression.FC$meanlogexprFC.shuffle <- NA
mean.expression.FC$medianlogexprFC.shuffle <- NA
mean.expression.FC$highlogexprFC.shuffle <- NA
mean.expression.FC$lowexprFC.shuffle <- NA
mean.expression.FC$meanexprFC.shuffle <- NA
mean.expression.FC$medianexprFC.shuffle <- NA
mean.expression.FC$highexprFC.shuffle <- NA
mean.expression.FC$propPosFC.shuffle <- NA
mean.expression.FC$propNegFC.shuffle <- NA

set.seed(1)
for (my.gene in unique(mean.expression.FC$gene)){
  for (my.cluster in unique(mean.expression.FC$cluster)){
    # For each gene cluster
    # I will use the mean expression evaluated from the pdf at each sample of the mcmc
    # To evaluate the mean expression and the fold-change
    print(paste(my.gene, my.cluster))
    mask <- mean.expression.FC$gene == my.gene &
      mean.expression.FC$cluster == my.cluster
    # I get the data
    means1 <- readLines(file.path(directory, mean.expression.FC$file[mask & mean.expression.FC$genotype != "wt"]))
    means2 <- readLines(file.path(directory, mean.expression.FC$file[mask & mean.expression.FC$genotype == "wt"]))
    # I put them in a dataframe where one of them was shuffle
    means <- data.frame(inv=as.numeric(means1), wt=sample(as.numeric(means2), length(means2), replace = F))
    # I evaluate the logFC and FC for each line = each sample pair
    means$logFC <- means$inv - means$wt
    means$FC <- delog(means$inv) / delog(means$wt)
    # I summarize them with mean/median/quantile
    mean.expression.FC$meanexpr.shuffle[mask & mean.expression.FC$genotype == "wt"] <- mean(delog(means$wt))
    mean.expression.FC$meanexpr.shuffle[mask & mean.expression.FC$genotype != "wt"] <- mean(delog(means$inv))
    mean.expression.FC$lowmeanexpr.shuffle[mask & mean.expression.FC$genotype == "wt"] <- quantile(delog(means$wt), probs = qm)
    mean.expression.FC$lowmeanexpr.shuffle[mask & mean.expression.FC$genotype!= "wt"] <- quantile(delog(means$inv), probs = qm)
    mean.expression.FC$highmeanexpr.shuffle[mask & mean.expression.FC$genotype == "wt"] <- quantile(delog(means$wt), probs = qp)
    mean.expression.FC$highmeanexpr.shuffle[mask & mean.expression.FC$genotype!= "wt"] <- quantile(delog(means$inv), probs = qp)
    mean.expression.FC$medianexpr.shuffle[mask & mean.expression.FC$genotype == "wt"] <- median(delog(means$wt))
    mean.expression.FC$medianexpr.shuffle[mask & mean.expression.FC$genotype!= "wt"] <- median(delog(means$inv))
    mean.expression.FC$meanlogexprFC.shuffle[mask] <- mean(means$logFC)
    mean.expression.FC$medianlogexprFC.shuffle[mask] <- median(means$logFC)
    mean.expression.FC$lowlogexprFC.shuffle[mask] <- quantile(means$logFC, probs = qm)
    mean.expression.FC$highlogexprFC.shuffle[mask] <- quantile(means$logFC, probs = qp)
    mean.expression.FC$meanexprFC.shuffle[mask] <- mean(means$FC)
    mean.expression.FC$lowexprFC.shuffle[mask] <- quantile(means$FC, probs = qm)
    mean.expression.FC$highexprFC.shuffle[mask] <- quantile(means$FC, probs = qp)
    mean.expression.FC$medianexprFC.shuffle[mask] <- median(means$FC)
    mean.expression.FC$propPosFC.shuffle[mask] <- sum(means$FC > 1) / nrow(means)
    mean.expression.FC$propNegFC.shuffle[mask] <- sum(means$FC < 1) / nrow(means)
  }
}
# I generate a table only with FC
simplified.expression.FC <- unique(mean.expression.FC[, c("cluster", "gene", grep("FC", colnames(mean.expression.FC), value = T))])
# I put the fold-change confidence interval to display it on the pdfs
simplified.expression.FC$value.label <- with(simplified.expression.FC, paste0(format(round(medianexprFC.shuffle, 2), nsmall = 2), "[", format(round(lowexprFC.shuffle, 2), nsmall = 2), "]",
                                                                              "^{", format(round(highexprFC.shuffle, 2), nsmall = 2), "}"))
# Figure 3Se
my.df <- subset(pdfs, file %in% log.evid.df$file[log.evid.df$is.best] & gene %in% c("Hoxd13", "Hoxa11", "Hoxd11", "Hoxa13") & cluster %in% paste0("cl", c(1, 6)))
my.df$cluster <- factor(my.df$cluster, levels = paste0("cl", c(1, 6)))

g <- ggplot(my.df, aes(x = x, color = genotype, fill = genotype)) +
  geom_path(aes(y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.3, color = NA) +
  facet_grid(cluster ~ gene, scales='free') +
  geom_text(data = simplified.expression.FC,
            inherit.aes = F,
            aes(label = value.label),
            x = 2.5, y = 2.5, parse = T) +
  theme_classic() +
  xlab("log(1+10'000*expression)") + 
  coord_cartesian(ylim=c(0, 3)) +
  expand_limits(x=-0.1) +
  ylab("Density") +
  ggtitle(my.model)
ggsave(paste0(directory.plot, "/fig3Se.pdf"), g, width = 8, height = 4.5)
