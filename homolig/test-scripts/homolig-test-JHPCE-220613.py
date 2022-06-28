import homolig as hg

inputDir = '/users/agirgis/homolig/homolig/data/test_data'
inputFile = inputDir + '/tcr_test.csv'

def main():
    adata = hg.homolig(input_file =  inputFile,
        seq_type = 'tcr',
        chains = 'paired',
        metric='aadist',
        species='human')

if __name__ == '__main__':
    main()

