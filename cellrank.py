import sys
import cellrank as cr
import scvelo as scv
import scanpy as sc
import warnings

sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2
warnings.simplefilter("ignore", category=UserWarning)

#load data from scvelo dynamic modeling processed data
adata = sc.read_h5ad('MAC_Plantir_pseudotime.h5ad')

pseudotime="palantir_pseudotime"

#run all lines
#getting started 
sc.pp.filter_genes(adata, min_cells=5)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.tl.pca(adata, random_state=0)
sc.pp.neighbors(adata, random_state=0)
sc.pl.embedding(adata, basis="umap", color=["clusters", pseudotime])
#Set up a kernel
from cellrank.kernels import PseudotimeKernel
pk = cr.kernels.PseudotimeKernel(adata, time_key=pseudotime)
#Let’s use this kernel to compute a cell-cell transition_matrix.
pk.compute_transition_matrix()

#Set up an estimator
from cellrank.estimators import GPCCA
g = GPCCA(pk)
adata.obs['clusters'] = adata.obs['clusters'].astype('category').values

#Identify initial & terminal states
g.fit(cluster_key="clusters")
g.predict_terminal_states(method="top_n", n_states=3)
g.predict_initial_states(allow_overlap=True)

#Compute fate probabilities and driver genes
g.compute_fate_probabilities()
g.plot_fate_probabilities(legend_loc="right")
cr.pl.circular_projection(adata, keys="clusters", legend_loc="right")

#Visualize expression trends
#separated by gene
#run calculation until Compute fate probabilities and driver genes
model = cr.models.GAMR(adata)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="MAGIC_imputed_data",
    genes=["Csf1r", "Axl", "Foxp1", "Ccr1", "Nampt", "Tgfbr2"], #MAC miR34a targets
    same_plot=True,
    ncols=3,
    time_key=pseudotime, #velocity_pseudotime or palantir_pseudotime or dpt_pseudotime
    hide_cells=True,
    legend_loc=None
)


