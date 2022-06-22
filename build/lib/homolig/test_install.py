import homolig as hig
print(dir(hig))
print(hig)

def main():
    '''
    adata = hg.homolig(input_file = 'bcr_test.csv', 
        seq_type = 'bcr', 
        chains = 'paired', 
        metric='aadist', 
        species='human')
    
    #adata = hg.homolig(input_file = 'datasets/vdjdb_formatted2.csv', 
    #adata = hg.homolig(input_file = 'homolig_paired_human.csv',
    '''
    #print(homolig.__version__)
    adata = hig.homolig(input_file = 'tcr_test.csv', 
        seq_type = 'tcr', 
        chains = 'paired',
        metric='aadist', 
        species='human')
    '''
    
    # WORKS
    adata = hg.homolig(input_file = 'test_seq.csv', 
        seq_type = 'seq', 
        chains = None, 
        metric='aadist', 
        species='human')
    
    adata = hg.homolig(input_file = 'CMV_M1_MART1_TSCAN_Donor_homolig.csv', 
        seq_type = 'tcr', 
        chains = 'beta', 
        metric='aadist', 
        species='human')
    
    adata = hg.homolig(input_file = 'references/emily_tcrdb_paired_human_imgt.csv', 
        seq_type = 'tcr', 
        chains = 'paired', 
        metric='aadist', 
        species='human')
        
    #print(adata.X)
    #print(adata.obs)
    '''
if __name__ == '__main__':
    main()
