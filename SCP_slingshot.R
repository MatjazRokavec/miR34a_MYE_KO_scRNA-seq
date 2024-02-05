library(SCP)
library(BiocParallel)
library(data.table)
register(SnowParam(workers = 14, progressbar = TRUE))

#load data
MAC <- readRDS(file = "MAC.rds")

#Trajectory inference (Slingshot)
MAC <- RunSlingshot(srt = MAC, group.by = "seurat_clusters", reduction = "UMAP", start="0")

