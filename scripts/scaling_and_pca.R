# scaling_and_pca.R
library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

pbmc <- readRDS(input_file)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Plots
pca_loadings <-VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave(filename = "results/pca_loading_plot.png", plot = pca_loadings, height = 6, width = 6, dpi = 300)

pca_dimplot <- DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave(filename = "results/pca_dimplot.png", plot = pca_dimplot, height = 6, width = 6, dpi = 300)

# Save single-dimension heatmap
png(filename = "results/pca_heatmap1.png", height = 1200, width = 1200, res = 300)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# Save multi-dimension heatmap
png(filename = "results/pca_heatmap15.png", height = 2000, width = 1200, res = 300)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Determine the ‘dimensionality’ of the dataset 
pca_elbowplot <- ElbowPlot(pbmc)
ggsave(filename = "results/elbow_plot.png", plot = pca_elbowplot, height = 6, width = 6, dpi = 300)

saveRDS(pbmc, file = output_file)
