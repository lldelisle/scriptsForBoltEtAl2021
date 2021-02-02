library(devtools)
install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("patchwork")
safelyLoadAPackageInCRANorBioconductor("ggplot2")

plot.folder <- "plots"
dir.create(plot.folder, showWarnings = F)
# We assume the working folder contains
# a folder for each sample
# the txt file meta_data_scRNAseq.txt (available in the same directory as this script)

# import individual cell-ranger output matrices.
folders <- grep("scRNA", list.dirs(), value = T)
names(folders) <- gsub("^./", "", folders)

genotypes <- unlist(lapply(strsplit(names(folders), "_"), "[[", 2))
names(genotypes) <- names(folders)

data <- lapply(folders, Read10X)
objects <- lapply(names(data), function(sample){
        # Create the seurat object
        seurat.obj <- CreateSeuratObject(counts = data[[sample]], project = sample, min.cells = 3, min.features = 500)
        # create filter for each object for mitochondrial genes (actual filtering done later)
        seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
        seurat.obj[["genotype"]] <- genotypes[sample]
        return(seurat.obj)
})
names(objects) <- names(data)
# check plot of features, counts, and percent mt reads
lapply(names(objects), function(sample){
        VlnPlot(objects[[sample]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        ggsave(file.path(plot.folder, paste0("QC_violins_", sample, ".png")), width = 7, height = 7)
        plot1 <- FeatureScatter(objects[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(objects[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        plot1 + plot2
        ggsave(file.path(plot.folder, paste0("QC_scatter_", sample, ".png")), width = 12, height = 7)
})

# subsetting only the cells with good qualities and reduction of doublets:
objects_sub <- lapply(objects, subset, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 5)

# The merge will not preserve reductions, graphs, logged commands, or feature-level metadata that were present in the original objects.
samples.merged <- merge(objects_sub[[1]],
                        subsetByNamesOrIndices(objects_sub, 2:length(objects_sub)),
                        add.cell.ids = names(objects_sub),
                        merge.data = TRUE,
                        project = "samples.merged")

table(samples.merged$orig.ident)
# E12.5_inv2_PPFL_scRNA_rep1 E12.5_inv2_PPFL_scRNA_rep2   E12.5_wt_PPFL_scRNA_rep1 
#                       3068                       3513                       3589 
table(samples.merged$genotype)
# inv2   wt 
# 6581 3589 

# We put wt first:
samples.merged$genotype <- factor(samples.merged$genotype, levels = c("wt", setdiff(unique(samples.merged$genotype), "wt")))

# normalization and scaling
samples.merged <- NormalizeData(samples.merged)
samples.merged <- FindVariableFeatures(samples.merged, selection.method = "vst")
all.genes <- rownames(samples.merged)


# Regress for cell-cycles genes
# details here https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html
# The updated genes list is cc.genes.updated.2019. Some warnings are normal. Ignore them. 
s.genes <- cc.genes.updated.2019$s.genes 
g2m.genes <- cc.genes.updated.2019$g2m.genes
samples.merged <- CellCycleScoring(samples.merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


head(samples.merged[[]])
# There should be a S.Score and G2M.Score with cell cycle phase columns added.
table(samples.merged$Phase)
#   G1  G2M    S 
# 5378 2519 2273 


# this will regress out the cell cycle genes.
samples.merged <- ScaleData(samples.merged, vars.to.regress = c("S.Score", "G2M.Score"), verbose = TRUE)


# check that the regression is good. can be done before and after cell cycle regression.
RidgePlot(samples.merged, features = c("Pcna", "Top2a", "Mcm6", "Mki67", "Hoxa13", "Hoxd13"), ncol = 2)
ggsave(file.path(plot.folder, "QC_cell_cycle.png"), width = 7, height = 7)

# this will run PCA for the whole data set
samples.merged <- RunPCA(samples.merged)


# Evaluate the PCs in different ways to identify the set to use. They recommend that you favor a larger number than smaller number of PCs.
png(file.path(plot.folder, "QC_elbow.png"))
ElbowPlot(samples.merged, ndims = 50, reduction = "pca")
dev.off()

# Other checks
# DimPlot(samples.merged, reduction = "pca", group.by = "Phase")
# DimPlot(samples.merged, reduction = "pca", group.by = "orig.ident")
# DimHeatmap(samples.merged, dims = 1:25, cells = 300, balanced = TRUE)


# UMAPing and clustering (considering number of dimensions from above)
samples.merged <- RunUMAP(samples.merged, dims = 1:25, seed.use = 22)
samples.merged <- FindNeighbors(samples.merged, dims = 1:25)
samples.merged <- FindClusters(samples.merged, resolution = 0.6)

DimPlot(samples.merged, reduction = "umap", split.by = "genotype", pt.size = 0.7, label = TRUE, label.size = 7, order = T)

# Because the clustering and the UMAP may change from machine, they are provided in a table
# Here are the commands used to generate the table:
# my.genes <- c("Hoxd13", "Hoxa11", "Hoxd11", "Hoxa13")
# df <- cbind(samples.merged[[]], 
#             samples.merged[["umap"]]@cell.embeddings, 
#             t(as.matrix(GetAssayData(object = samples.merged, assay = "RNA", slot = "counts")[my.genes,])))
# write.table(df, "meta_data_scRNAseq.txt", quote = F, sep = "\t")

# To apply the data in the table to your seurat object:
meta.data <- read.delim("meta_data_scRNAseq.txt", row.names = 1)[, c("seurat_clusters", "UMAP_1", "UMAP_2")]
# Adjust the clusters
samples.merged$seurat_clusters <- factor(meta.data[colnames(samples.merged), "seurat_clusters"],
                                         levels = sort(unique(samples.merged$seurat_clusters)))
Idents(samples.merged) <- samples.merged$seurat_clusters

# We put the UMAP coordinates
samples.merged[["umap"]]@cell.embeddings <- as.matrix(meta.data[colnames(samples.merged), c("UMAP_1", "UMAP_2")])

DimPlot(samples.merged, reduction = "umap", split.by = "genotype", pt.size = 0.7, label = TRUE, label.size = 7, order = T) & xlim(c(-8, 8)) & ylim(c(-6, 6))

# To remove the clusters that are outside the main group
mainCluster <- subset(samples.merged, idents = c("10", "12", "13", "14", "15", "16"), invert = TRUE)

# Figure 3c: UMAP
p1 <- DimPlot(mainCluster, reduction = "umap", split.by = "genotype", pt.size = 0.7, label = T, label.size = 7, order = T, ncol = 2, cols = "Paired")
p1 + coord_cartesian(xlim = c(-8, 8), ylim = c(-6, 6))
ggsave(file.path(plot.folder, "Figure3c.pdf"), width = 12, height = 7)

my.fixed.colors <- c("magenta4", "mediumblue", "darkgreen", "orange", "magenta")
names(my.fixed.colors) <- c("Hoxd13", "Hoxd11", "Hoxa11", "Hoxa13", "Shox2")
# Figure 3d: 
my.genes <- c("Hoxd13", "Hoxa11")
for (i in 1:length(my.genes)){
        FeaturePlot(samples.merged, features = my.genes[i], order = T, pt.size = 1.2, label = T,
                    label.size = 7, ncol = 1, cols = c("lightgrey", my.fixed.colors[my.genes[i]]), split.by = "genotype",
                    max.cutoff = 2.5) + 
                theme(legend.position = 'right') &
                coord_cartesian(xlim = c(-8, 8), ylim = c(-6, 6)) 
        ggsave(file.path(plot.folder, paste0("Figure3d_", my.genes[i], ".pdf")), width = 10, height = 5)
}

# Figure 3e: Violin plots in cluster 1 and 6
my.clusters.nb <- c(1, 6)
my.clusters.colors <- c("dodgerblue4", "orange")
for (i in 1:length(my.clusters.nb)){
        my.cluster <- subset(samples.merged, subset = seurat_clusters == my.clusters.nb[i])
        VlnPlot(my.cluster, features = c("Hoxd13", "Hoxa11", "Hoxd11", "Hoxa13"), group.by = "genotype",
                y.max = 3, flip = T, pt.size = 0.4,
                cols = c("grey", my.clusters.colors[i]), ncol = 4)
        ggsave(file.path(plot.folder, paste0("Figure3e_cluster", my.clusters.nb[i], ".pdf")), width = 10, height = 3)
}

# Figure S3b: selected genes in each cluster
VlnPlot(mainCluster, features = c("Hoxd13", "Hoxa11","Hoxd11", "Hoxa13", "Shox2"), 
        group.by = "seurat_clusters", split.by = "genotype", same.y.lims = T, stack = T, flip = T) +
        scale_fill_manual(values=c('wt'="grey", 'inv2'="pink"))
ggsave(file.path(plot.folder, "Figure3Sb.pdf"), width = 12, height = 7)


wt_subset <- subset(samples.merged, subset = genotype == "wt")
inv2_subset <- subset(samples.merged, subset = genotype == "inv2")

# Figure S3c: selected genes on UMAP
my.genes <- c("Hoxd13", "Hoxd11", "Hoxa11", "Hoxa13", "Shox2")
for (i in 1:length(my.genes)){
        p2 <- FeaturePlot(wt_subset, features = my.genes[i], order = T, pt.size = 1.2, label = T,
                          label.size = 7, ncol = 1, cols = c("lightgrey", my.fixed.colors[my.genes[i]])) +
                ggtitle('wt')
        p3 <- FeaturePlot(inv2_subset, features = my.genes[i], order = T, pt.size = 1.2, label = T,
                          label.size = 7, ncol = 1, cols = c("lightgrey",  my.fixed.colors[my.genes[i]])) + 
                ggtitle('inv2')
        g <- (p2 + coord_cartesian(xlim = c(-8, 8), ylim = c(-6, 6)) | p3 + coord_cartesian(xlim = c(-8, 8), ylim = c(-6, 6))) + plot_layout(guides = 'collect')
        ggsave(file.path(plot.folder, paste0("Figure3Sc_", my.genes[i], ".pdf")), g, width = 10, height = 5)
}

# Figure 3Sd: proportion of cells with detected expression 
my.genes <- c("Hoxd13", "Hoxa11", "Hoxd11", "Hoxa13")
fixed.colors.a11.a13.d11 <- c(my.fixed.colors, "grey")
names(fixed.colors.a11.a13.d11) <- c(names(my.fixed.colors), "No expression")
big.table <- read.delim("meta_data_scRNAseq.txt", row.names = 1)
big.table$Hoxd13_status <- "NoHoxd13"
big.table$Hoxd13_status[big.table$Hoxd13 > 0] <- "Hoxd13"
my.genesx <- setdiff(my.genes, c("Hoxd13"))
for (genex in my.genesx){
        big.table[, paste0(genex, "_status")] <- "No expression"
        big.table[, paste0(genex, "_status")][big.table[, genex] > 0] <- genex
}

tables <- do.call(rbind, lapply(my.genesx, function(genex){
        df <- as.data.frame(table(big.table$seurat_clusters, big.table[, paste0(genex, "_status")], big.table$genotype, big.table$Hoxd13_status))
        df$genex <- genex
        return(df)
}))
colnames(tables)[1:4] <- c("cluster", "category", "genotype", "Hoxd13_status")
tables$cluster <- as.character(tables$cluster)
# I want to add 2 special clusters
temp.table <- ddply(subset(tables, cluster %in% c(1, 6)),
                    .(category, genotype, genex, Hoxd13_status),
                    summarise,
                    Freq = sum(Freq))
temp.table$cluster <- "1,6"
tables <- rbind(tables, temp.table[, colnames(tables)])
temp.table <- ddply(subset(tables, cluster %in% c(3,4,5)),
                    .(category, genotype, genex, Hoxd13_status),
                    summarise,
                    Freq = sum(Freq))
temp.table$cluster <- "3,4,5"
tables <- rbind(tables, temp.table[, colnames(tables)])

# Evaluate proportion
tables <- ddply(tables, .(genex, cluster, genotype, Hoxd13_status), mutate,
                prop = Freq / sum(Freq))
# Evaluate the pval
tables <- ddply(tables, .(genex, cluster, genotype), mutate,
                fisherpval = fisher.test(matrix(c(Freq[category != "No expression" & Hoxd13_status == "Hoxd13"],
                                                  Freq[category != "No expression" & Hoxd13_status == "NoHoxd13"],
                                                  Freq[category == "No expression" & Hoxd13_status == "Hoxd13"],
                                                  Freq[category == "No expression" & Hoxd13_status == "NoHoxd13"]),
                                                nrow = 2))$p.value)
tables$category <- factor(tables$category, levels = names(fixed.colors.a11.a13.d11))
tables$Hoxd13_status <- factor(tables$Hoxd13_status, levels = c("NoHoxd13", "Hoxd13"))

for (my.genex in my.genesx){
        my.tables <- subset(tables, genex == my.genex & ((cluster %in% c("3", "4", "5", "3,4,5", "1", "6", "1,6"))))
        my.tables$cluster <- factor(my.tables$cluster, levels = c("3", "4", "5", "3,4,5", "1", "6", "1,6"))
        ggplot(my.tables, aes(x = Hoxd13_status, y=prop)) +
                geom_bar(aes(fill=category), stat = 'identity') +
                geom_text(data = unique(subset(my.tables, select = c("cluster", "genex", "fisherpval", "genotype"))),
                          aes(label = format(fisherpval, digits = 2)),
                          x = 2, y = 0.5, size = 2.5) +
                facet_grid(genotype ~ cluster, labeller = label_both) +
                scale_fill_manual(values=fixed.colors.a11.a13.d11) +
                ggtitle(paste0("Correlation ", my.genex, " Hoxd13"))
        ggsave(file.path(plot.folder,paste0("Figure3Sd_", my.genex, ".pdf")),
               height = 5, width = 12)
}
