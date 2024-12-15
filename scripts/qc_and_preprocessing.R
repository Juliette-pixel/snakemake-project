# qc_and_preprocessing.R
library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

pbmc <- readRDS(input_file)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
QC_violin <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "results/QC_violin_plot.png", plot = QC_violin, height = 6, width = 12, dpi = 300)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
feature_scatter <- plot1 + plot2
ggsave(filename = "results/feature_scatter_plot.png", plot = feature_scatter, height = 6, width = 12, dpi = 300)


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
saveRDS(pbmc, file = output_file)


