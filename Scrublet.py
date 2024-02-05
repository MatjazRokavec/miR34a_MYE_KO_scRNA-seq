%matplotlib inline
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys
import scanpy as sc
#import read10x

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

counts_matrix = scipy.io.mmread('LYS_CRE_NEW/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes('LYS_CRE_NEW/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))
counts_matrix = scipy.io.mmread('wt_NEW/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes('wt_NEW/filtered_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

dbl_rate = counts_matrix.shape[0]/1000 * 0.008
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbl_rate, sim_doublet_ratio = 2)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=3, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

scrub.plot_histogram();
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True);

results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
dataframe = pd.concat([results, scores], axis=1)
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

dataframe.to_csv(os.path.join('scrublet_results.tsv'), sep = "\t", index = False)

