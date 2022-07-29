input_file=F:/homolig/homolig/data/test_data/bcr-paired-test.csv
seq_type=bcr
chains=paired
metric=aadist
species=human

python runHomolig.py -i $input_file -s $seq_type -c $chains -m $metric -p $species
