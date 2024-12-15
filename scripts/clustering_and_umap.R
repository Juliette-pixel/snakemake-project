# clustering_and_umap.R
library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

pbmc <- readRDS(input_file)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

umap_dimplot <- DimPlot(pbmc, reduction = "umap")
ggsave(filename = "results/umap_dimplot.png", plot = umap_dimplot, height = 6, width = 6, dpi = 300)


saveRDS(pbmc, file = output_file)
