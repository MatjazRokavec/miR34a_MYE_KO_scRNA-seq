library(CellChat)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

combi_wt <- readRDS(file = "combi_wt.rds")

#create a cellchat object from seurat object
DefaultAssay(combi_wt) <- "RNA" #RNA data
cellchat <- createCellChat(object = combi_wt, group.by = "ident")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel; doesnt work
cellchat <- identifyOverExpressedGenes(cellchat) #takes long
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat) #takes long
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat, file = "cellchat_wt.rds")

combi_CRE <- readRDS(file = "combi_CRE.rds")
#create a cellchat object from seurat object
DefaultAssay(combi_CRE) <- "RNA" #RNA data
cellchat <- createCellChat(object = combi_CRE, group.by = "ident")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel; doesnt work
cellchat <- identifyOverExpressedGenes(cellchat) #takes long
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat) #takes long
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
saveRDS(cellchat, file = "cellchat_CRE.rds")

cellchat_wt <- readRDS(file = "cellchat_wt.rds")
cellchat_CRE <- readRDS(file = "cellchat_CRE.rds")
object.list <- list(wt = cellchat_wt, CRE = cellchat_CRE)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), merge.data = TRUE)
saveRDS(cellchat, file = "cellchat_combi.rds")




#show only specific pathways
pathways.show <- c("CSF")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway
gg1 <- netAnalysis_contribution(cellchat_wt, signaling = pathways.show)
gg2 <- netAnalysis_contribution(cellchat_CRE, signaling = pathways.show)
gg1+gg2
