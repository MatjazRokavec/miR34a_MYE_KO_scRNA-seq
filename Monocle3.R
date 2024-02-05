#requires packages leidenbase>0.1.25, igraph>1.5.0
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(data.table)

MAC <- readRDS(file = "MAC.rds")
cds <- as.cell_data_set(MAC)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
cds <- cluster_cells(cds, resolution=0.0005)
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)

#Learn the trajectory graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#order cells according to pseudotime
cds <- order_cells(cds) #manually choose starting cell
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
