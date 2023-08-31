#/usr/bin/env python

import os
import sys
import homolig as homolig
from optparse import OptionParser
import time

def main(options):
    input_file=options.input_file
    seq_type=options.seq_type
    chains=options.chains
    metric=options.metric
    species=options.species
    homolig_out = homolig.homolig(input_file, 
                                  seq_type, 
                                  chains, 
                                  metric, 
                                  species)
    
if __name__ == '__main__':
    
    parser = OptionParser()

    parser.add_option("-i","--input_file", dest="input_file",
                      help="the input file")
    parser.add_option("-s", "--seq_type", dest="seq_type",
                      help="the seq type ['tcr','bcr', or 'seq']")
    parser.add_option("-c", "--chains", dest="chains",
                      help="chains ['alpha', 'beta', 'paired', 'light', 'heavy', or None]")
    parser.add_option("-m", "--metric", dest="metric",
                      help="metrics ['aadist', 'aadist_euc', or 'blosum62']",default="aadist"),
    parser.add_option("-p", "--species", dest="species",
                      help="species ['human', 'mouse', 'mas-night-monkey', 'rhesus-monkey','alpaca', 'bovine', 'camel', 'catfish', 'chicken', 'chondrichthyes''cod', 'crab-eating-macaque', 'dog', 'dolphin', 'ferret','goat', 'gorilla', 'horse', 'naked-mole-rat','nonhuman-primates', 'pig', 'platypus', 'rabbit', 'rat','salmon', 'sheep', 'teleostei', 'trout', 'zebrafish']",default="human")
    
    (options, args) = parser.parse_args()

    main(options)
