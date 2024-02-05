#Change environment to:
#Change environmnet, Open Tools -> preferences -> Python interpreter -> Use the following Python interpreter:
#D:/Programs/python38/python.exe

#pip install -U scvelo

import scanpy as sc
import scvelo as scv
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


adata = sc.read_h5ad('MAC.h5ad')

adata
# plot umap to check
sc.pl.umap(adata, color=['clusters'], frameon=False, legend_loc='on data', title='')

# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('wt_NEW.loom', cache=True)
ldata2 = scv.read('LYS_CRE_NEW.loom', cache=True)

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
ldata2.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()

# concatenate the three loom
ldata = ldata1.concatenate([ldata2])

#clean
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color='clusters', frameon=False, legend_loc='on data', title='')

#preprocess data
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=20000) #original = 2000
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=20000) #warnings
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#velocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata) #takes long time; 

#project velocity
scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters')
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120, color='clusters')

#proportion of spliced and unspliced
scv.pl.proportions(adata, groupby='clusters') #works after you run preprocess data and run #velocity

#Speed and coherence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95]) 

df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

#Velocity graph and pseudotime
scv.pl.velocity_graph(adata, threshold=.1)
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax) 
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot') 

#Dotplot, cluster specific genes
sc.tl.rank_genes_groups(adata,'clusters', use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=10, groupby='clusters', use_raw=False, dendrogram=False)


adata.write('MAC_scvelo.h5ad', compression='gzip')

