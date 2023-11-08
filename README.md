<!-- ABOUT THE PROJECT -->
## What is Homolig?

Homolig is a python package to physiochemically compare immune receptor and
epitope sequences.

<!-- GETTING STARTED -->
## Getting Started

### Installation

To get a local copy of Homolig up and running follow these simple example steps.
(Once publicly available, this will be installable directly from pip without the
need for cloning.) 

1. Clone the repository
```bash
git clone https://github.com/FertigLab/homolig
```
2. Pip install
```bash
cd homolig 
python3 -m pip install .
```

You should now be able to import homolig within python. Alternatively, you may call a wrapper script from the terminal as shown below (recommended). 

## Documentation
Homolig uses the IMGT/V-QUEST reference directory release 202214-2 (05 April
2022).

## File format
The input file is a comma separated file containing the TCR or BCR CDR3 amino acid sequence and varible
gene name in IMGT format.   

For cases where paired alpha and beta chain information is available:  
| CDR3.beta.aa | TRBV | CDR3.alpha.aa | TRAV |
| --- | --- |  --- | --- |
| CASSAGTSPTDTQYF | TRBV6-4*01 | CAVMDSSYKLIF | TRAV1-2*01 |

### Basic Usage
Recommended usage is through the wrapper homolig_wrapper.py, located in ./homolig/. 
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

To cluster the results of a Homolig run, you may use the wrapper clusterHomolig.py: 

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

To generate a UMAP reduction based on the NxN distance matrix, you may use the wrapper homolig_write_umap.py: 

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

## Citation

If you use this software, please cite our manuscript:

<!-- CONTACT -->
## Contact

Please send feedback to Alex Girgis -
<agirgis3@jhu.edu>


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.
