import scanpy as sc
import pandas as pd
import numpy as np
import diopy
import anndata as ad
import scvi
import os
import argparse
import configparser
# parse cmdline arguments
parser = argparse.ArgumentParser(description='get arguments')
parser.add_argument('-c', '--config', dest = "configuration", help='configuration file with absolute path')
args = parser.parse_args()
configuration = args.configuration

sc.logging.print_header()
config = configparser.ConfigParser()
config.read(configuration)
config.sections()
scvi_out = config["Path"]['scvi_out']
merge_out = config["Path"]['merge_out']
sample_info = config['Path']['sample_info']
merge_result = config['Result_file']['merge_result']
batch_key = config['Integration']['batch_key']
ntop = int(config['Integration']['ntop'])
nlayer = int(config['Integration']['nlayer'])
nlatent = int(config['Integration']['nlatent'])

print('check1')

if not os.path.exists(os.path.join(scvi_out)):
    os.makedirs(os.path.join(scvi_out))

adata = sc.read_h5ad(os.path.join(merge_out, merge_result))
adata.obs_names_make_unique()

adata.layers['hvgcounts']=adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw=adata
sc.pp.highly_variable_genes(
	adata
	, flavor='seurat'
	, n_top_genes=ntop
	, subset=True
	, layer='hvgcounts'
	, batch_key=batch_key
)

print('check2')

scvi.data.setup_anndata(adata, layer='hvgcounts', batch_key=batch_key)
vae=scvi.model.SCVI(adata, n_layers=nlayer, n_latent=nlatent)
vae.train()
vae.save(scvi_out + 'scvi_model')
adata.obsm['X_scVI']=vae.get_latent_representation()

print('check3')

sc.pp.neighbors(adata, use_rep='X_scVI', random_state=1234)
sc.tl.leiden(adata, resolution=1, random_state=1234)
sc.tl.umap(adata, random_state=1234)

adata.write(os.path.join(scvi_out, 'scvi.h5ad'))

bdata = adata.raw.to_adata()
# sc.pp.neighbors(bdata, use_rep='X_scVI', random_state=1234)
# sc.tl.leiden(bdata, resolution=1, random_state=1234)
# sc.tl.umap(bdata, random_state=1234)
bdata.write('/storage/chenlab/Users/qingnanl/humanre/out/scvi_out/scvi.raw.h5ad')
diopy.output.write_h5(adata = bdata,
                      file = '/storage/chenlab/Users/qingnanl/humanre/out/scvi_out/scvi.raw.h5')




