#Run k-means on Homolig similarity matrices.
#Using mbkmeans on raw similarity matrices - not doing any kind of dimensional reduction.
#User inputs 'expected' number of clusters C.  This script will sweep across that value 1/2C, C, 2C, 4C
#AAG 27 May 2023

import sys
import argparse
from time import localtime, strftime
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans

#ARGPARSE------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', help = "Input .h5ad file")
parser.add_argument('-c', '--num_clusters', help = "Expected number of clusters. May be any integer.")
parser.add_argument('-o', '--output',  help = "Desired output file path/filename. Defaults to input file directory.")
parser.add_argument('-r', '--report',  help = "Whether to generate silhouette score report. Logical.")
#Read arguments from command line.
args = parser.parse_args()

if args.input is None:
    sys.exit('Please specify input file.')
if args.num_clusters is None:
    sys.exit('Please specify expected number of groups in data.')
if args.output is None:
    args.output = args.input[0:-5] + '_clusters.csv'
if args.report is None:
    args.report = False

#Load in the data /assign variables. ---------------------------------------------

inputfile=args.input
outputfile=args.output
C=int(args.num_clusters)
r = args.report
now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('[' + now + ']', 'Homolig version 0.1 mbkmeans Clustering')  # Find a way to add __version__ attribute to package at later date.
print('                         input_file: ', inputfile)
print('                        output_file: ', outputfile)
print('                       num_clusters: ', C)
print('             make silhouette report: ', r)


#Execute clustering!--------------------------------------------------------------

adata = sc.read_h5ad(inputfile)

x = np.asarray(adata.X.todense())


kmeans1 = MiniBatchKMeans(n_clusters = round(C / 2), random_state = 99, batch_size = 1024, n_init = 10)
kmeans2 = MiniBatchKMeans(n_clusters = C, random_state = 99, batch_size = 1024, n_init = 10)
kmeans3 = MiniBatchKMeans(n_clusters = 2 * C, random_state = 99, batch_size = 1024, n_init = 10)
kmeans4 = MiniBatchKMeans(n_clusters = 4 * C, random_state = 99, batch_size = 1024, n_init = 10)


clusters1 = kmeans1.fit_predict(x)
clusters2 = kmeans2.fit_predict(x)
clusters3 = kmeans3.fit_predict(x)
clusters4 = kmeans4.fit_predict(x)

info_list = list(zip( clusters1, clusters2, clusters3, clusters4))

df = pd.DataFrame( info_list,
                   columns = [ 'C1', 'C2', 'C3', 'C4'])

obs = adata.obs
obs.reset_index(drop=True, inplace=True)
df.reset_index(drop=True,inplace=True)
cluster_info = pd.concat([obs, df], axis=1)

cluster_info.to_csv(outputfile)
#If indicated, draft report plots of silhouette score.---------------
if args.report:
    print('Report is True.')
    #Make name for cluster output file.
    report_name = args.input[0:-5] + '_cluster-silhouette-report.png'
    #Use example code to draft graphs of silhouette score.
    #Consider a basic dim red to graph.

    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)

        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

        # Initialize the clusterer with n_clusters value and a random generator
        # seed of 10 for reproducibility.
        clusterer = KMeans(n_clusters=n_clusters, n_init="auto", random_state=10)
        cluster_labels = clusterer.fit_predict(X)

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(X, cluster_labels)
        print(
            "For n_clusters =",
            n_clusters,
            "The average silhouette_score is :",
            silhouette_avg,
        )

        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(X, cluster_labels)

        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.nipy_spectral(float(i) / n_clusters)
            ax1.fill_betweenx(
                np.arange(y_lower, y_upper),
                0,
                ith_cluster_silhouette_values,
                facecolor=color,
                edgecolor=color,
                alpha=0.7,
            )

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])



#Print Closing Message and Exit-------------------------------------
now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('[' + now + ']', 'Homolig clustering completed.')
