import sys
import argparse
import homolig as hg
import warnings
warnings.filterwarnings("ignore")
#24 September 2022 AAG
#Updated 2 January 2023

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', help = "Input .csv file")
parser.add_argument('-s', '--seq', help = "Sequence type. May be one of [tcr, bcr, seq].")
parser.add_argument('-c', '--chains', help = "Chain locus. May be one of [alpha, beta, heavy, light]. Can be omitted if --seq == 'seq'. ")
parser.add_argument('-m', '--metric', help = "Distance matrix used in comparisons. Default is aadist.")
parser.add_argument('-sp', '--species', help = "Species for which to query V gene sequences. Default is human.")
parser.add_argument('-mode', '--mode', help = "Either 'pairwise' sequence comparison or 'axb' between two sequence groups.")
parser.add_argument('-i2', '--input2', help = "Second sequence group with which to compare first file.")
parser.add_argument('-o', '--output',  help = "Desired output file path/filename. Defaults to input file directory.")
parser.add_argument('-v', '--verbose', help= "Specify verbosity during execution.")
parser.add_argument('-g', '--save_germline', help= 'Save CDR alignments separately during execution.')
parser.add_argument('-si', '--save_reformatted_input', help= 'Save input file after renaming V genes (may be useful in post-analysis).')
#Read arguments from command line.
args = parser.parse_args()

#print(args)

if args.input is None:
    sys.exit('Please specify input file.')
if args.seq is None:
    sys.exit('Please specify --seq as one of [tcr, bcr, seq]')
if args.chains is None:
    if (args.seq != "'seq'"):
        sys.exit('Please specify --chains as one of [alpha, beta, heavy, light, paired]')
    args.chains = 'beta'
if args.metric is None:
    args.metric = 'aadist'
if args.species is None:
    args.species = 'human'
if args.mode is None:
    args.mode = 'pairwise'
if args.save_germline is None:
    args.save_germline = False
args.verbose = True
if args.save_reformatted_input is None: 
    args.save_reformatted_input = False
#No problem if args.input2 or output is None, so will not check.


args.input = args.input.replace("'", "")
args.seq = args.seq.replace("'", "")
args.chains = args.chains.replace("'", "")
args.metric = args.metric.replace("'", "")
args.species = args.species.replace("'", "")
args.mode = args.mode.replace("'","")
if args.input2 is not None:
    args.input2 = args.input2.replace("'","")



#Now that we have finished processing inputs, let's execute the main homolig function.

def main():
    adata = hg.homolig_format_checker(input_file= args.input,
        seq_type= args.seq,
        chains = args.chains,
        metric= args.metric,
        species= args.species,
        mode= args.mode,
        input2= args.input2,
        output_file= args.output,
        verbose= args.verbose, 
        save_reformatted_input = args.save_reformatted_input
                       )

if __name__ == '__main__':
    main()
