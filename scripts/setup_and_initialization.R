# setup_and_initialization.R
library(dplyr)
library(sp)
library(SeuratObject)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
output_file <- args[2]

pbmc.data <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
saveRDS(pbmc, file = output_file)
