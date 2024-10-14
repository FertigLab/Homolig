
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from scipy.sparse import csr_matrix
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import umap

#3 October 2022

#ARGPARSE========
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', help = "Path to .h5ad file containing cluster classifications.")
parser.add_argument('-o', '--output', help = "Path to write output csv file.")

args = parser.parse_args()

args.input = args.input.replace("'", "")
args.output = args.output.replace("'", "")

#Read inputs and calculate UMAP============
adata = sc.read_h5ad(args.input)

nComp = min(100, adata.shape[0]) #not totally positive if shape[0] is the correct dimension here.
svd = TruncatedSVD(n_components=nComp, n_iter=7, random_state=42)
svd.fit(adata.X)
print('SVD Axis Variances: ', svd.explained_variance_ratio_)
print('Total Variance Explained: ', round(sum(svd.explained_variance_ratio_),3))

svd_dat = svd.transform(adata.X)

#Now get a UMAP.
umapObj = umap.UMAP()
umap_dat = umapObj.fit_transform(svd_dat)

#Now write those coordinates to file.==========
udf = pd.DataFrame(umap_dat)
udf.index = adata.obs_names
udf = udf.rename(columns={0: 'UMAP1', 1: 'UMAP2'})
udf.to_csv(args.output)
