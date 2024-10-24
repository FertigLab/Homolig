#Run k-means on Homolig similarity matrices.
#Using mbkmeans on raw similarity matrices - not doing any kind of dimensional reduction.
#User inputs 'expected' number of clusters C.  This script will sweep across that value 1/2C, C, 2C, 4C
#AAG 27 May 2023
#Added silhouette score computation 11 September 2023 AAG
import sys
import argparse
from time import localtime, strftime
import warnings
warnings.filterwarnings("ignore")

from anndata import read_h5ad
import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_samples
from sklearn.decomposition import PCA
#ARGPARSE------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', help = "Input .h5ad file")
parser.add_argument('-c', '--num_clusters', help = "Expected number of clusters. May be any integer.")
parser.add_argument('-o', '--output',  help = "Desired output file path/filename. Defaults to input file directory.")
parser.add_argument('-d', '--num_components',  help = "Number of principal components kept for clustering. Defaults to 100.")

#Read arguments from command line.
args = parser.parse_args()

if args.input is None:
    sys.exit('Please specify input file.')
if args.num_clusters is None:
    sys.exit('Please specify expected number of groups in data.')
if args.output is None:
    args.output = args.input[0:-5] + '_' + str(args.num_clusters) + '_clusters.csv'
if args.num_components is None:
    args.num_components = 100

#Load in the data /assign variables. ---------------------------------------------

inputfile=args.input
outputfile=args.output
C=int(args.num_clusters)
nPC=int(args.num_components)

now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('[' + now + ']', 'Homolig version 0.4 mbkmeans Clustering')  # Find a way to add __version__ attribute to package at later date.
print('                         input_file: ', inputfile)
print('                        output_file: ', outputfile)
print('                       num_clusters: ', C)
print('                     num_components: ', nPC)


#Execute clustering!--------------------------------------------------------------

adata = read_h5ad(inputfile)

x = np.asarray(adata.X.todense())
pca = PCA(n_components=nPC)
pca.fit(x)
P = pca.transform(x)
print('           Explained Variance Ratio: ', round(sum(pca.explained_variance_ratio_), 4))
kmeans1 = MiniBatchKMeans(n_clusters = C, random_state = 99, batch_size = 1024, n_init = 10)

clusters1 = kmeans1.fit_predict(P)


info_list = list(zip( clusters1))

df = pd.DataFrame( info_list,
                   columns = [ 'Cluster'])
#Silhouette score computation
df['sil_score'] = silhouette_samples(x, df['Cluster'])


#Prepare and save final output
obs = adata.obs
obs.reset_index(drop=True, inplace=True)
df.reset_index(drop=True,inplace=True)
cluster_info = pd.concat([obs, df], axis=1)

cluster_info.to_csv(outputfile)


#Print Closing Message and Exit-------------------------------------
now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('[' + now + ']', 'Homolig clustering completed.')
