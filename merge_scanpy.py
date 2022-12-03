import scanpy as sc
import pandas as pd
import numpy as np
import diopy
import anndata as ad
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
pdir = config['Path']['pdir']
out_dir = config['Path']['out_dir']
QC_out = config['Path']['QC_out']
merge_out = config["Path"]['merge_out']
sample_info = config['Path']['sample_info']
merge_result = config['Result_file']['merge_result']

if not os.path.exists(os.path.join(out_dir)):
    os.makedirs(os.path.join(out_dir))

# if not os.path.exists(os.path.join(pdir, out_dir, QC_out)):
#    os.makedirs(os.path.join(pdir, out_dir, QC_out))

if not os.path.exists(os.path.join(merge_out)):
    os.makedirs(os.path.join(merge_out))

smpl = pd.read_csv(sample_info)# make sure the csv file is comma-delimited
print("sample file")
print(smpl)
adatas = []

for index, row in smpl.iterrows():
  dir = QC_out + row['unique'] + '.h5'
  print(dir)
  adata = diopy.input.read_h5(file = dir)
  adata.var_names_make_unique()
  adatas.append(adata)

#    
merged_adata = ad.concat(adatas, join='outer')
merged_adata.write(os.path.join(merge_out, merge_result))