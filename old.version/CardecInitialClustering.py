# please refer to cardec tutorial https://github.com/jlakkis/CarDEC
import pandas as pd
import os
import numpy as np
import pickle
from copy import deepcopy
from shutil import move
import warnings

"""Machine learning and single cell packages"""
import sklearn.metrics as metrics
import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
"""CarDEC Package"""
from CarDEC import CarDEC_API


# run CarDEC
results_file = '~/human_retina.cardec.15_donor.h5ad'
adata = sc.read("~/human_retina_update.h5ad")
#test with subsampled; then go with real data 
#sc.pp.subsample(adata, n_obs=10000, random_state=0, copy=False)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
CarDEC = CarDEC_API(adata,
        batch_key = "donor", n_high_var = 3000, LVG = True,normalize_samples = False, log_normalize = False)

CarDEC.build_model(n_clusters = 15)
CarDEC.make_inference()
CarDEC.model_counts()
print(CarDEC.dataset)
CarDEC.dataset.write(results_file)

cardec = sc.read('~/human_retina.cardec.15_donor.h5ad')
denoised = deepcopy(cardec.layers['denoised counts'][:, :])
denoised_round=np.round(denoised, 2)
new = AnnData(denoised_round)
new.obs=cardec.obs
new.var=cardec.var
new.uns=cardec.uns
new.obsm=cardec.obsm
q = deepcopy(new.obsm['cluster memberships']) #The cluster membership numpy array
labels = np.argmax(q, axis=1)
labels = [str(x) for x in labels]
new.obs["cluster"] = list(labels)

sc.pp.neighbors(new, n_neighbors = 15, use_rep = 'embedding')
sc.tl.umap(new)
#
sc.pl.umap(new, color=['cluster'])
plt.savefig('./cardec.15_donor.cardec.cluster.tiff', dpi = 300)
#sc.pl.umap(new, color=['donor','group'])
#plt.savefig('./cardec.15_donor.cardec.dimred.tiff')
sc.pl.umap(new, color=['donor'], legend_fontsize = '5')
plt.savefig('./cardec.15_donor.by.donor.tiff', dpi = 300)

# checking a few markers
sc.pl.umap(new, color=["PDE6A","ARR3", "RLBP1", "GFAP"])
plt.savefig('./cardec.15_donor.cardec.marker1.tiff', dpi = 300)
sc.pl.umap(new, color=["GAD2","SLC6A9", "RBPMS", "C1QB"])
plt.savefig('./cardec.15_donor.cardec.marker2.tiff', dpi = 300)
sc.pl.umap(new, color=[ "GRIK1", "GRM6","ONECUT2", "RPE65"])
plt.savefig('./cardec.15_donor.cardec.marker3.tiff', dpi = 300)

#some clusters contain more than one cell types; try leiden to separate clusters
sc.tl.leiden(new, resolution = 0.3)
sc.pl.umap(new, add_outline=True,legend_loc='on data',color=['leiden'])
plt.savefig('./cardec.15_donor.cardec.leiden.tiff')
sc.tl.rank_genes_groups(new, 'leiden', method='t-test')

sc.pl.rank_genes_groups(new, n_genes=25, sharey=False)
plt.savefig('./cardec.15.marker.leiden.tiff')

#update leiden


new.obs['annotate_major'] = (
    new.obs["leiden"]
    .map(lambda x: {'0':'Rod', '1':'Rod', '2':'NN', '3':'AC','4':'Cone', '5':'BC',
        '6':'HC', '7':'RGC', '8':'AC', '9':'BC', '10':'BC', '11':'BC',
        '12':'AC', '13':'AC','14':'AC', '15':'AC','16':'BC', '17':'BC',
        '18':'AC', '19':'AC','20':'AC', '21':'BC','22':'AC', '23':'BC',
        '24':'BC', '25':'RGC','26':'BC', '27':'Rod','28':'NN', '29':'NN'}.get(x, x))
    .astype("category")
)

results_file = './CarDEC.annotated.h5ad'
new.write(results_file)


adata = sc.read('/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/human_retina_update.h5ad')
adata.obs['annotate_major']=new.obs['annotate_major']
#split and output
AC = adata[adata.obs["annotate_major"] == "AC"]
result_AC="/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/cardec_output/cluster_15_donor/AC/AC.h5ad"
AC.write(result_AC)

BC = adata[adata.obs["annotate_major"] == "BC"]
result_BC="/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/cardec_output/cluster_15_donor/BC/BC.h5ad"
BC.write(result_BC)
#....
