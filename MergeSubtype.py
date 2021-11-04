import pandas as pd
import os
import numpy as np
import pickle
from copy import deepcopy
from shutil import move
import warnings

import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt

Rod = sc.read("./Rod/Rod_annotated.h5ad")
Rod = Rod.raw.to_adata()
Rod.obs["major_type"] = "Rod"
Rod.obs["ml_class"] = "Rod"
Rod.obs["subtype"] = "Rod"
Rod.write("./Rod/Rod_annotated.h5ad")


BC = sc.read("./BC_13/BC_annotated.h5ad")
BC.obs["major_type"] = "BC"
BC.obs['ml_class'] = (
    BC.obs["cluster"]
    .map(lambda x: {'0':'OFFBC', '1':'ONBC', '2':'ONBC', '3':'OFFBC','4':'ONBC',
    '5':'OFFBC', '6':'ONBC', '7':'OFFBC', '8':'OFFBC','9':'ONBC',
    '10':'ONBC', '11':'ONBC', '12':'OFFBC'}.get(x, x))
    .astype("category")
)

BC.obs['subtype'] = (
    BC.obs["cluster"]
    .map(lambda x: {'0':'BC0', '1':'BC1', '2':'BC2', '3':'BC3','4':'BC4',
    '5':'BC5', '6':'BC6', '7':'BC7', '8':'BC8','9':'BC9',
    '10':'BC10', '11':'BC11', '12':'BC12'}.get(x, x))
    .astype("category")
)
BC.write("./BC_13/BC_annotated.h5ad")


....

data = Rod.concatenate(Cone, AC, BC, HC,
     RGC, MG, NN, join = 'outer')

data.var['gene_ids'] = data.var.index
data.var['feature_types'] = "Gene Expression"
# data.var has too many columns due to outer join; remove them
df = data.var[data.var.columns[354:356]]
data.var = df
data.write("./data_annotated.h5ad")
