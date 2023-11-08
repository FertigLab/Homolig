
import argparse
import pandas as pd
import scanpy as sc
import numpy as np
import warnings
warnings.filterwarnings("ignore")
#Take h5ad output from homolig and write matrix as csv along with supporting obs/var metadata.
#AAG 2 January 2023

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = "Input .h5ad file")
parser.add_argument('-o', '--output',  help = "Desired output file path/filename. Defaults to input file directory.")

#Read arguments from command line.
args = parser.parse_args()
args.input = args.input.replace("'", "")
if args.output is not None:
    args.output = args.input.replace("'", "")
else:
    args.output = args.input[:-5] + '.csv'

infile = args.input
outfile = args.output
metafile = outfile[:-4] + '_sequence_info.csv'
adata = sc.read_h5ad(infile)

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

input_agg.to_csv(metafile, index=False)
mat.to_csv(outfile)
