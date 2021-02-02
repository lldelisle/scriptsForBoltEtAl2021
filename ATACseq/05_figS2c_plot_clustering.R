options(stringsAsFactors = F)

# Install dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("reshape")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggdendro")
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# One of the dependency of ggpubr is only available for R=4.
if ( ! require("ggpubr") & ! rversionAbove(4)){
  packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  safelyLoadAPackageInCRANorBioconductor("ggpubr")
}

setwd(commandArgs(TRUE)[1])

ggfin <- list()
for (position in c("left", "right")){
  data <- read.delim(paste0("Figure2Sd_", position, "_data.txt"), quote = "\'")
  data[is.na(data)] <- 0
  # Combine the three first columns
  my.data <- cbind(apply(data[, 1:3], 1, function(v){paste0(v[1], ":", v[2], "-", v[3])}), data[, 4:ncol(data)])
  # Remove the extension
  # colnames(my.data) <- gsub("_ATAC.bw", "", colnames(my.data))
  colnames(my.data) <- sapply(colnames(my.data), function(s){paste(strsplit(s, "_")[[1]][3:2], collapse = "\n")})
  colnames(my.data)[1] <- "region"
  # Reverse the order
  my.data$region <- factor(my.data$region, levels = sort(unique(my.data$region), decreasing = T))
  # melt
  melt.df <- melt(my.data)
  # Do a pearson correlation
  cor.m <- cor(my.data[, 2:ncol(my.data)])
  # Clustering
  my.clust <- hclust(as.dist(1 - cor.m))
  dd <- as.dendrogram(my.clust)
  my.clust2 <- as.hclust(reorder(dd, 1:length(my.clust$labels), agglo.FUN = mean))
  g1 <- ggdendrogram(my.clust2) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_blank(),
          text = element_text(size = 15),
          plot.margin = margin(0, 10, 0, 40, unit = "pt"))
  # Now plot the data
  my.gg <- list()
  for (v in my.clust2$labels[my.clust2$order]){
    if (grepl("DFL", v)){
      my.colors <- c("#fff5f0", "#67000d")
    }
    if (grepl("PFL", v)){
      my.colors <- c("#f7fbff", "#08306b")
    }
    if (grepl("FB", v)){
      my.colors <- c("#ffffff", "#000000")
    }
    my.gg[[v]] <- ggplot(subset(melt.df, variable == v), aes(x = variable, y = region)) +
      geom_tile(aes(fill = value)) +
      scale_fill_gradient(low = my.colors[1], high = my.colors[2], 
                          n.breaks = 4, limits = c(0, NA)) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = 'bottom',
            legend.title = element_blank(),
            axis.text.x = element_text(size = 20),
            legend.text = element_text(size = 15),
            plot.margin = margin(0, 0, 0, 0, unit = "pt")) +
      ylab("") +
      xlab("")
  }
  g2 <- ggarrange(plotlist = my.gg, nrow = 1)
  ggfin[[position]] <- ggarrange(g1, g2, nrow = 2, heights = c(1, 4))
}
x <- 4
ggexport(ggarrange(plotlist = ggfin), filename = "Figure2Sd.pdf", width = 3 * x, height = 1.5 * x)
