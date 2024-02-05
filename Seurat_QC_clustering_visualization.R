library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(EnhancedVolcano)


# Load the 10x datasets (feature matrix)
wt.data <- Read10X(data.dir = "wt/filtered_feature_bc_matrix")
LYS_CRE.data <- Read10X(data.dir = "LYS_CRE/filtered_feature_bc_matrix")
wt <- CreateSeuratObject(counts = wt.data, project = "wt")
LYS_CRE <- CreateSeuratObject(counts = LYS_CRE.data, project = "LYS_CRE")

#remove cells with less than 200 nFeature_RNA and more than 25% mitochondrial RNA
wt$percent.MT <- PercentageFeatureSet(wt,pattern="^mt-")
LYS_CRE$percent.MT <- PercentageFeatureSet(LYS_CRE,pattern="^mt-")
VlnPlot(wt, features=c("nCount_RNA", "nFeature_RNA", "percent.MT"))
VlnPlot(LYS_CRE, features=c("nCount_RNA", "nFeature_RNA", "percent.MT"))
wt <- subset(wt, nFeature_RNA > 200 & percent.MT < 25)
LYS_CRE <- subset(LYS_CRE, nFeature_RNA > 200 & percent.MT < 25)

#remove cell coublets found with Scrublet
doublets <- read.delim(file = "doublets_scrublet_wt.txt", sep="\t", dec=".")
doublets <- doublets$wt
wt <- wt[,!colnames(wt) %in% doublets]
doublets <- read.delim(file = "doublets_scrublet_CRE.txt", sep="\t", dec=".")
doublets <- doublets$CRE
LYS_CRE <- LYS_CRE[,!colnames(LYS_CRE) %in% doublets]

#generate a list of wt and LYS_CRE objects
combined <- list(wt=wt,LYS_CRE=LYS_CRE)

# normalize and identify variable features for each dataset independently
combined <- lapply(X = combined, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000) #LogNormalize (default) or CLR: centered log ratio normalization (margin = 1 (per gene) or 2 (per cell)); RC: relative counts (for counts per million set scale.factor = 1e6)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = combined)
#Perform Integration
immune.anchors <- FindIntegrationAnchors(object.list = combined, anchor.features = features)
to_integrate <- Reduce(intersect, lapply(immune.anchors@object.list, rownames))
combi <- IntegrateData(anchorset = immune.anchors, features.to.integrate = to_integrate)

#Clustering
DefaultAssay(combi) <- "integrated"
all.genes <- rownames(combi)
combi <- ScaleData(combi, verbose = FALSE, features=all.genes)
combi <- RunPCA(combi, features=VariableFeatures(object=combi), verbose = FALSE)
combi <- RunUMAP(combi, reduction = "pca", dims = 1:50)
combi <- FindNeighbors(combi, reduction = "pca", dims = 1:50)
combi <- FindClusters(combi, resolution = 1) #0.5 - 1.5
saveRDS(combi, file = "combi.rds")

#Plots
DefaultAssay(combi) <- "RNA"
DimPlot(combi, reduction = "umap", split.by = "orig.ident", label = T) #group by clusters and make two plots (wt and LYS_CRE)
markers.to.plot <- c("Epcam", "Cdh1", "Ocln", #epi
                     "Adgre1", "Cd68", "Csf1r", "Lyz2", #MAC; Adgre1 = F4/80
                     "Csf3r", "Trem1", "Cxcr2",  #Neutrophil; also Retnlg and Lcn2 are good neutro markers
                     "Tpsab1", "Cma1", "Tph1", #Mast cells
                     "Cd69", "Il3ra", "Cpa3", #Basophil
                     "Clec9a", "Batf3", #DC
                     "Cd3d", "Cd3e", #T-cells
                     "Cd79a", "Cd79b", #B-cells
                     "Col6a2", "Col3a1", #Fibro
                     "Plvap", "Vwf", #Endo
                     "Hbb-bt", "Hba-a1") #erythro
DotPlot(combi, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() + coord_flip()
VlnPlot(combi, features = markers.to.plot, stack=TRUE, flip=TRUE, pt.size = 0, split.by = "orig.ident", split.plot=TRUE) #separated by clusters

# find markers for every cluster compared to all remaining cells
clusters <- combi
DefaultAssay(clusters) <- "RNA"
clusters.markers <- FindAllMarkers(clusters, only.pos = TRUE)
clusters.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
fwrite(clusters.markers, file = "markers.txt", sep="\t", dec=".", row.names=TRUE)
#Heatmap with top 20 markers for each cluster
clusters.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 20) %>%
  ungroup() -> top10
DoHeatmap(Tcell, features = top10$gene) + NoLegend() + scale_fill_gradientn(colours = c("lightblue","white","red"))
        
#Differentially expressed genes
DefaultAssay(combi) <- "RNA"
combi$celltype.CRE <- paste(Idents(combi), combi$orig.ident, sep="_")
combi$celltype <- Idents(combi)
Idents(combi) <- "celltype.CRE"
markers <- FindMarkers(combi, ident.1 = 'MAC_wt', ident.2 = 'MAC_LYS_CRE', logfc.threshold = 0)
fwrite(markers, file = "DEG.txt", sep="\t", dec=".", row.names=TRUE)

#volcano plot
markers <- read.delim(file = "volcano.txt", sep="\t", dec=".", row.names = 1)
keyvals.col <- ifelse(markers$miR34a.score < 500, 'blue',
                      ifelse(markers$published > 0, 'red', ''))
names(keyvals.col)[keyvals.col == 'blue'] <- 'predicted miR-34a target genes'
names(keyvals.col)[keyvals.col == 'red'] <- 'published miR-34a target genes'
p1 <- EnhancedVolcano(markers,
                      lab = rownames(markers),
                      x = 'rev_log2FC',
                      y = 'p_val',
                      colCustom = keyvals.col,
                      colAlpha = 1, #transparency
                      labSize = 0,
                      pointSize = 2,
                      vline = c(0))
p1 + ggplot2::xlab("avg. log2 fold change Lys_CRE/wt")

#export data to create a h5ad file
seurat_obj <- combi
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
seurat_obj$Cell_type <- seurat_obj@active.ident
write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='data')
writeMM(counts_matrix, file='counts.mtx')
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

