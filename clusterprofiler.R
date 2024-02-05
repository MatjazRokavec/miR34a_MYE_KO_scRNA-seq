library(clusterProfiler)
library(enrichplot)
library(msigdbr)


# Read in the data; create a list of a txt file; text file consists of different groups with Entrez IDs for each gene
x <- scan("mac.txt", what="", sep="\n")

y <- strsplit(x, "[[:space:]]+")
names(y) <- sapply(y, `[[`, 1)
y <- lapply(y, `[`, -1)

#compare cluster to GO, WP, KEGG
#ck <- compareCluster(a, fun = enrichKEGG, organism="hsa")
ck <- compareCluster(geneCluster = y, fun = enrichGO, OrgDb='org.Mm.eg.db',
                     pvalueCutoff = 0.000000001,pAdjustMethod = "BH", qvalueCutoff = 0.000000001)
ck <- compareCluster(geneCluster = y, fun = enrichWP, organism="mmu")

#enrich cluster in any Msigdb gene set
#load Msigdb gene sets
msigdbr_show_species()
#H: hallmark gene sets
#C1: positional gene sets
#C2: curated gene sets
#C3: motif gene sets
#C4: computational gene sets
#C5: GO gene sets
#C6: oncogenic signatures
#C7: immunologic signatures
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% #Hallmark
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG") %>% #KEGG
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% #Cell type
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>% #Oncogenic gene sets
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% #Immunologic gene sets
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% #GO_Biological+Process
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC") %>% #GO_Cellular_Component
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF") %>% #GO_Molecular_Function
  dplyr::select(gs_name, entrez_gene)
ck <- compareCluster(geneCluster = y, fun = enricher, TERM2GENE=m_t2g,
                     pvalueCutoff = 0.05,pAdjustMethod = "BH", qvalueCutoff = 0.05)
ck <- compareCluster(geneCluster = y, fun = enricher, TERM2GENE=m_t2g,
                     pvalueCutoff = 0.0001,pAdjustMethod = "BH")
ck <- compareCluster(geneCluster = y, fun = enricher, TERM2GENE=m_t2g,
                     pvalueCutoff = 0.000000001,pAdjustMethod = "BH")
ck <- compareCluster(geneCluster = y, fun = enricher, TERM2GENE=m_t2g,
                     pvalueCutoff = 0.0000000000000000001,pAdjustMethod = "BH", qvalueCutoff = 0.0000000000000000001)
ck <- compareCluster(geneCluster = y, fun = enricher, TERM2GENE=m_t2g,
                     pvalueCutoff = 0.000000000000000000000000000001,pAdjustMethod = "BH", qvalueCutoff = 0.00000000000000000000000000001)

ck@compareClusterResult$qvalue <- -log10(ck@compareClusterResult$qvalue)

#plot
dotplot(ck, showCategory=10, color="qvalue") + scale_color_gradient(low="sienna1", high="red2")

