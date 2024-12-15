# differential_markers_and_visualization.R
library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)



args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
markers_output <- args[2]
plot_output <- args[3]

pbmc <- readRDS(input_file)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

marker_violin <- VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
ggsave("results/marker_violin.png", marker_violin, height = 6, width = 6, dpi = 300)

marker_violin_raw <- VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE)
ggsave("results/marker_violin_raw.png", marker_violin_raw, height = 6, width = 6, dpi = 300)

marker_feature <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
ggsave("results/marker_feature.png", marker_feature, height = 8, width = 10, dpi = 300)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
marker_heatmap <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("results/marker_heatmap.png", marker_heatmap, height = 8, width = 10, dpi = 300)

saveRDS(pbmc.markers, file = markers_output)

plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18),
        legend.text = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = plot_output, height = 7, width = 12, plot = plot, quality = 50)
