<!-- ABOUT THE PROJECT -->
## What is Homolig?

Homolig is a python package to physiochemically compare immune receptor and
epitope sequences. There are two modules:
(1) Compute pairwise sequence distance between sequences / Assign Clusters / Write UMAP
(2) Generate descriptive statistics on physiochemical properties of selected sequences. 

<!-- GETTING STARTED -->
## Getting Started

### Installation

The below steps are necessary to use the Homolig clustering module, written in Python3. 
1. Clone the repository
```bash
git clone https://github.com/FertigLab/Homolig
```
2. If using docker: may build local docker image using Dockerfile. Alternatively, install Homolig within a virtual environment. A global installation is not recommended due to several package version dependencies.

```bash
#apt install python3-venv  #if not already installed
cd Homolig 
python3 -m venv ./homolig-venv
source ./homolig-venv/bin/activate
python3 -m pip install .
```

3. Now try using Homolig to cluster example immune repertoires. If successful, output files from this test script will populate the repository parent directory. 

```bash
# sudo chmod +x ./tests/testcode-python-clustering-module.sh #if you don't have privileges to execute test script
./tests/testcode-python-clustering-module.sh
```


## Documentation
Homolig uses the IMGT/V-QUEST reference directory release 202214-2 (05 April
2022).

## File format
The input file is a comma separated file containing the TCR or BCR CDR3 amino acid sequence and varible
gene name in IMGT format. See example inputs in ./test-data/.  

For cases where paired alpha and beta chain information is available:  
| CDR3.beta.aa | TRBV | CDR3.alpha.aa | TRAV |
| --- | --- |  --- | --- |
| CASSAGTSPTDTQYF | TRBV6-4*01 | CAVMDSSYKLIF | TRAV1-2*01 |

### Basic Usage: Pairwise Distances
Recommended usage for Module 1: Pairwise similarity scores is through the wrapper homolig_wrapper.py, located in ./homolig/. 
```python
python3 $WRAPPERDIR/homolig_wrapper.py  -h
usage: homolig_wrapper.py [-h] [-i INPUT] [-s SEQ] [-c CHAINS] [-m METRIC] [-sp SPECIES] [-mode MODE] [-i2 INPUT2] [-o OUTPUT] [-v VERBOSE]
                          [-g SAVE_GERMLINE] [-si SAVE_REFORMATTED_INPUT]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input .csv file
  -s SEQ, --seq SEQ     Sequence type. May be one of [tcr, bcr, seq].
  -c CHAINS, --chains CHAINS
                        Chain locus. May be one of [alpha, beta, heavy, light]. Can be omitted if --seq == 'seq'.
  -m METRIC, --metric METRIC
                        Distance matrix used in comparisons. Default is aadist.
  -sp SPECIES, --species SPECIES
                        Species for which to query V gene sequences. Default is human.
  -mode MODE, --mode MODE
                        Either 'pairwise' sequence comparison or 'axb' between two sequence groups.
  -i2 INPUT2, --input2 INPUT2
                        Second sequence group with which to compare first file.
  -o OUTPUT, --output OUTPUT
                        Desired output file path/filename. Defaults to input file directory.
  -v VERBOSE, --verbose VERBOSE
                        Specify verbosity during execution.
  -g SAVE_GERMLINE, --save_germline SAVE_GERMLINE
                        Save CDR alignments separately during execution.
  -si SAVE_REFORMATTED_INPUT, --save_reformatted_input SAVE_REFORMATTED_INPUT
                        Save input file after renaming V genes (may be useful in post-analysis).

```

To cluster the results of a Homolig run, you may use the wrapper clusterHomolig.py. In the instance of large datasets, you may consider clusterHomolig.pca.py which generates a UMAP using a PCA reduction of the output similarity scores: 

```python
python3 $WRAPPERDIR/clusterHomolig.py  -h
usage: clusterHomolig.py [-h] [-i INPUT] [-c NUM_CLUSTERS] [-o OUTPUT]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input .h5ad file
  -c NUM_CLUSTERS, --num_clusters NUM_CLUSTERS
                        Expected number of clusters. May be any integer.
  -o OUTPUT, --output OUTPUT
                        Desired output file path/filename. Defaults to input file directory.

```

To generate a UMAP reduction based on the NxN similarity matrix, you may use the wrapper homolig_write_umap.py: 

```python
python3 $WRAPPERDIR/clusterHomolig.py  -h
usage: clusterHomolig.py [-h] [-i INPUT] [-c NUM_CLUSTERS] [-o OUTPUT]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input .h5ad file
  -c NUM_CLUSTERS, --num_clusters NUM_CLUSTERS
                        Expected number of clusters. May be any integer.
  -o OUTPUT, --output OUTPUT
                        Desired output file path/filename. Defaults to input file directory.
```
### Basic Usage: Descriptive Module 
Functions to describe the physiochemical properties of arbitrary sequence groups are written in R. 
Please see functions in ./homolig/_rcode/score-sequences.r. For a walkthrough of basic usage see ./tests/testcode-r-characterization-module.rmd. 

## Citation

Upon publication, please cite our manuscript, "Comparative Assessment of Physiochemical Metrics for the Clustering of Adaptive Immune Receptor Repertoires" (currently in preparation). For now, please cite this repository. 

Girgis, A., Davis-Marcisak, E., & Palmer, T. Homolig - tools for the physiochemical analysis of immune repertoires [Computer software]

<!-- CONTACT -->
## Contact

Please send feedback to Alex Girgis -
<agirgis3@jhu.edu>


<!-- LICENSE -->
## License

Distributed under AGPL 3.0. See `LICENSE` for more information.
