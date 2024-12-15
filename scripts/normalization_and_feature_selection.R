# normalization_and_feature_selection.R
library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

pbmc <- readRDS(input_file)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
variable_features <- plot1 + plot2
ggsave(filename = "results/variable_features_plot.png", plot = variable_features, height = 6, width = 12, dpi = 300)

saveRDS(pbmc, file = output_file)


