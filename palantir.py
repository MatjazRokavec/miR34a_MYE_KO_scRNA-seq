import palantir
import scanpy as sc
import pandas as pd
import os

# Plotting
import matplotlib
import matplotlib.pyplot as plt
# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)
# Inline plotting
%matplotlib inline

#load data
adata = sc.read_h5ad('MAC.h5ad')

# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(adata)

# Use scanpy functions to visualize umaps or FDL
sc.pl.umap(adata, color=["clusters"])

#Diffusion maps visualization
palantir.plot.plot_diffusion_components(adata)
ms_data = palantir.utils.determine_multiscale_space(adata)

#MAGIC imputation
imputed_X = palantir.utils.run_magic_imputation(adata)
#plot MAGIC imputed data
sc.pl.embedding(
    adata,
    basis="umap",
    layer="MAGIC_imputed_data",
    color=["Csf1r", "Axl", "Foxp1", "Ccr1", "Nampt", "Tgfbr2"], #MAC
    frameon=False,
    ncols=2,

)

#compute pseudotime
start_cell = "AACACACGTCATCAGT-1_1" #barcode of an early cell
pr_res = palantir.core.run_palantir(
    adata, start_cell, num_waypoints=500, use_early_cell_as_start = True
)
palantir.plot.plot_palantir_results(adata, s=3)

adata.write('MAC_palantir_pseudotime.h5ad', compression='gzip')












