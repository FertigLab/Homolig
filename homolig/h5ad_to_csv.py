
import argparse
from time import localtime, strftime
import pandas as pd
from anndata import read_h5ad
import numpy as np
import warnings
warnings.filterwarnings("ignore")
#Take h5ad output from homolig and write matrix as csv along with supporting obs/var metadata.
#AAG 2 January 2023

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = "Input .h5ad file")
parser.add_argument('-o', '--output',  help = "Desired output file path/filename. Defaults to input file directory.")
parser.add_argument('-m', '--meta',  help = "Whether to save metadata file (h5ad obs). Defaults to False.")
#Read arguments from command line.
args = parser.parse_args()
args.input = args.input.replace("'", "")
if args.output is not None:
    args.output = args.input.replace("'", "")
else:
    args.output = args.input + '.csv'
if args.meta is None:
    args.meta = False

infile = args.input
outfile = args.output
saveMeta = args.meta
metafile = outfile[:-4] + '_sequence_info.csv'


now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('[' + now + ']', 'Homolig version 0.2 Write Matrix to CSV')  # Find a way to add __version__ attribute to package at later date.
print('                         input_file: ', infile)
print('                        output_file: ', outfile)
print('                          save_meta: ', saveMeta)


adata = read_h5ad(infile)

if adata.var.shape[1] > 0: #If we used axb mode to generate input, and input contains both var and obs pd DFs
    mat = pd.DataFrame(adata.X)
    mat.columns = adata.var['Homolig.ID']
    mat.index = adata.obs['Homolig.ID']
    input_agg = pd.concat([adata.obs, adata.var])
else: #If we used pairwise mode to generate input, and input only contains obs DF.
    mat = pd.DataFrame(adata.X.todense())
    mat.columns = adata.obs['Homolig.ID']
    mat.index = adata.obs['Homolig.ID']
    input_agg = adata.obs

if metafile is True:
	input_agg.to_csv(metafile, index=False)

mat = np.round(mat,3)	
mat.to_csv(outfile, header=False,index=False)


#Print Closing Message and Exit-------------------------------------
now = strftime("%Y-%m-%d %H:%M:%S", localtime())
print('[' + now + ']', 'Homolig Write Matrix to CSV completed.')
