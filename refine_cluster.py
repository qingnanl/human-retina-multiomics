import igraph as ig
import leidenalg as la
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors, kneighbors_graph
from sklearn.metrics import pairwise_distances, silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.model_selection import ParameterGrid, train_test_split, cross_val_score
import os
import umap
# import umap.plot
# import matplotlib.pyplot as plt
# import seaborn as sns
from sklearn import svm
import argparse

# python refine_cluster.py -p /storage/chenlab/Users/qingnanl/humanre/out/subtypes/ \
#                          -i RGC_subset_emb.csv \
#                          -o RGC_subset_df.csv \
#                          -c RGC_subset_cluster.csv\
# parse cmdline arguments
parser = argparse.ArgumentParser(description='get arguments')

parser.add_argument('-p', '--pdir', dest = "pdir", help='parent directory for input')
parser.add_argument('-i', '--input', dest = "inputdata", help='input dataframe in csv format')
parser.add_argument('-o', '--output', dest = "output", help='output dataframe')
parser.add_argument('-c', '--cluster', dest = "cluster", help='cluster information')
args = parser.parse_args()

pdir = args.pdir
inputdata = args.inputdata
output = args.output
cluster_info = args.cluster

emb = pd.read_csv(pdir + inputdata, index_col=0)
n_cells = len(emb.index)

def find_clusters(embeddings, resoscan=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                  nnscan=[10, 20, 30, 40, 50], min_f1 = 0.95):
      pwdt = pairwise_distances(embeddings)
      np.fill_diagonal(pwdt, 0)
      clf = svm.SVC(kernel='linear')
      param_grid = {'resoscan': resoscan, 'nnscan': nnscan}
      dict_list = list(ParameterGrid(param_grid))
      n_search = len(dict_list)
      df = pd.DataFrame(columns=["nn", "resolution", "silhouette", "CH", "DB", "f1", "n_cluster"])
      for i in range(n_search):
            p = dict_list[i].get("nnscan")
            resolution = dict_list[i].get("resoscan")
            A = kneighbors_graph(embeddings.to_numpy(), p, mode='connectivity', include_self=True)
            g = ig.Graph.Adjacency((A.toarray() > 0).tolist())
            # Add edge weights and node labels.
            g.es['weight'] = A[A.nonzero()]
            g.vs['label'] = embeddings.index
            partition = la.find_partition(g, la.RBConfigurationVertexPartition,
                                          resolution_parameter=resolution, seed=1)
            labels = partition.membership
            avg_silhouette = silhouette_score(pwdt, labels, metric="precomputed")
            CH = calinski_harabasz_score(embeddings, labels)
            DB = davies_bouldin_score(embeddings, labels)
            scores = cross_val_score(clf, embeddings, np.asarray(labels), cv=3, scoring='f1_macro')
            n_cluster = len(set(labels))
            values_to_add = {"nn": p, "resolution": resolution, "silhouette": avg_silhouette,
                             "CH": CH, "DB": DB, "f1": scores.mean(), "n_cluster": n_cluster}
            row_to_add = pd.Series(values_to_add, name=i)
            df = df.append(row_to_add)

      index_combination = df['silhouette'].idxmax()
      p1 = df.iloc[index_combination]['nn']
      resolution1 = df.iloc[index_combination]['resolution']
      A = kneighbors_graph(embeddings.to_numpy(), int(p1), mode='connectivity', include_self=True)
      g = ig.Graph.Adjacency((A.toarray() > 0).tolist())
      g.es['weight'] = A[A.nonzero()]
      g.vs['label'] = embeddings.index
      partition = la.find_partition(g, la.RBConfigurationVertexPartition,
                                    resolution_parameter=resolution1, seed=1)
      labels = partition.membership
      
      return p1, labels, df
      
nn, labels, out_df = find_clusters(embeddings=emb)
umap_embedding = umap.UMAP(random_state=1, n_neighbors=int(nn), metric='correlation').fit_transform(emb)
umap_df = pd.DataFrame(umap_embedding, columns = ['UMAP1','UMAP2'])
umap_df['label'] = labels
umap_df.to_csv(pdir + cluster_info)
out_df.to_csv(pdir + output)







