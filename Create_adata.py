"""Use the cardecenv conda environment"""
"""This is to read 10x files, combine them and get their source info (donor/region/group) added"""
import pandas as pd
import os
import numpy as np
import warnings
import scanpy as sc
from anndata import AnnData

results_file = '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/human_retina_update.h5ad'
D19D014_macular = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D014_macular/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D014_macular.obs['group'] = 'D19D014_macular'
D19D014_macular.obs['donor'] = 'D19D014'
D19D014_macular.obs['region'] = 'macular'
D19D014_macular.obs['treatment'] = 'macular'
#
D19D013_macular = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D013_mac/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D013_macular.obs['group'] = 'D19D013_macular'
D19D013_macular.obs['donor'] = 'D19D013'
D19D013_macular.obs['region'] = 'macular'
D19D013_macular.obs['treatment'] = 'macular'
#
D19D015_macular = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D015_mac/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D015_macular.obs['group'] = 'D19D015_macular'
D19D015_macular.obs['donor'] = 'D19D015'
D19D015_macular.obs['region'] = 'macular'
D19D015_macular.obs['treatment'] = 'macular'
#
D19D016_macular = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D016_mac/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D016_macular.obs['group'] = 'D19D016_macular'
D19D016_macular.obs['donor'] = 'D19D016'
D19D016_macular.obs['region'] = 'macular'
D19D016_macular.obs['treatment'] = 'macular'
#
D19D014_fovea = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D014_fovea/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D014_fovea.obs['group'] = 'D19D014_fovea'
D19D014_fovea.obs['donor'] = 'D19D014'
D19D014_fovea.obs['region'] = 'fovea'
D19D014_fovea.obs['treatment'] = 'fovea'
#
D19D014_foveaR = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D014_foveaR/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D014_foveaR.obs['group'] = 'D19D014_fovea'
D19D014_foveaR.obs['donor'] = 'D19D014'
D19D014_foveaR.obs['region'] = 'fovea'
D19D014_foveaR.obs['treatment'] = 'fovea'
#
D19D013_fovea = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D013_fovea/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D013_fovea.obs['group'] = 'D19D013_fovea'
D19D013_fovea.obs['donor'] = 'D19D013'
D19D013_fovea.obs['region'] = 'fovea'
D19D013_fovea.obs['treatment'] = 'fovea'
#
D19D013_foveaR = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D013_foveaR/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D013_foveaR.obs['group'] = 'D19D013_fovea'
D19D013_foveaR.obs['donor'] = 'D19D013'
D19D013_foveaR.obs['region'] = 'fovea'
D19D013_foveaR.obs['treatment'] = 'fovea'
#
#
D19D015_fovea = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D015_fovea/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D015_fovea.obs['group'] = 'D19D015_fovea'
D19D015_fovea.obs['donor'] = 'D19D015'
D19D015_fovea.obs['region'] = 'fovea'
D19D015_fovea.obs['treatment'] = 'fovea'
#
D00012_NeuNM = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/00012_NeuNM/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D00012_NeuNM.obs['group'] = 'D00012_per'
D00012_NeuNM.obs['donor'] = 'D00012'
D00012_NeuNM.obs['region'] = 'per'
D00012_NeuNM.obs['treatment'] = 'NeuNM'
#
D00012_NeuNT = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/00012_NeuNT/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D00012_NeuNT.obs['group'] = 'D00012_per'
D00012_NeuNT.obs['donor'] = 'D00012'
D00012_NeuNT.obs['region'] = 'per'
D00012_NeuNT.obs['treatment'] = 'NeuNT'

#
D17D13_NeuNT = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/17D13_NeuNT/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D17D13_NeuNT.obs['group'] = 'D17D13_NeuNT'
D17D13_NeuNT.obs['donor'] = 'D17D13'
D17D13_NeuNT.obs['region'] = 'per'
D17D13_NeuNT.obs['treatment'] = 'NeuNT'
#
D19D013_NeuNT = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D013_NeuNT/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D013_NeuNT.obs['group'] = 'D19D013_per'
D19D013_NeuNT.obs['donor'] = 'D19D013'
D19D013_NeuNT.obs['region'] = 'per'
D19D013_NeuNT.obs['treatment'] = 'NeuNT'
#
D19D014_NeuNT = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D014_NeuNT/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D014_NeuNT.obs['group'] = 'D19D014_per'
D19D014_NeuNT.obs['donor'] = 'D19D014'
D19D014_NeuNT.obs['region'] = 'per'
D19D014_NeuNT.obs['treatment'] = 'NeuNT'

#
D19D015_NeuNT = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D015_NeuNT/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D015_NeuNT.obs['group'] = 'D19D015_per'
D19D015_NeuNT.obs['donor'] = 'D19D015'
D19D015_NeuNT.obs['region'] = 'per'
D19D015_NeuNT.obs['treatment'] = 'NeuNT'
#
D19D015_NeuNM = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D015_NeuNM/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D015_NeuNM.obs['group'] = 'D19D015_per'
D19D015_NeuNM.obs['donor'] = 'D19D015'
D19D015_NeuNM.obs['region'] = 'per'
D19D015_NeuNM.obs['treatment'] = 'NeuNM'

#
D19D015_per = sc.read_10x_mtx(
    '/storage/singlecell/qingnanl/human_10x_process/37_snRNAseq_re/data_output_folder/post_QC_data/19D015_per/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
D19D015_per.obs['group'] = 'D19D015_per'
D19D015_per.obs['donor'] = 'D19D015'
D19D015_per.obs['region'] = 'per'
D19D015_per.obs['treatment'] = 'per'



adata = D19D014_macular.concatenate(D19D013_macular, D19D015_macular, D19D016_macular, D19D014_fovea,
     D19D014_foveaR, D19D013_fovea, D19D013_foveaR, D19D015_fovea, D00012_NeuNM, D00012_NeuNT,  D17D13_NeuNT,
     D19D013_NeuNT, D19D014_NeuNT, D19D015_NeuNM, D19D015_NeuNT, D19D015_per, join = 'outer') # please note 'outer' here.
adata.write(results_file)
