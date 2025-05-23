import sys
import os
import pandas as pd
import numpy as np
import more_itertools as mit # get sliding windows
from Bio import SeqIO
from Bio.Seq import translate
import itertools
import anndata as ad
from tqdm.contrib.concurrent import process_map
from multiprocessing import Pool
from multiprocessing import cpu_count
from scipy import sparse
from homolig import substitution_matrices
import homoligcpp as homoligcpp
from time import localtime, strftime
import pkg_resources
DATA_PATH = pkg_resources.resource_filename('homolig', 'data/')

def _species_lookup(species):
    species_dict = {'mas-night-monkey':'Aotus_nancymaae',
                 'rhesus-monkey':'Macaca_mulatta',
                 'alpaca':'Vicugna_pacos',
                 'bovine':'Bos_taurus',
                 'camel':'Camelus_dromedarius',
                 'cat':'Felis_catus',
                 'catfish':'Ictalurus_punctatus',
                 'chicken':'Gallus_gallus',
                 'chondrichthyes':'chondrichthyes',
                 'cod':'Gadus_morhua',
                 'crab-eating-macaque':'Macaca_fascicularis',
                 'dog':'Canis_lupus_familiaris',
                 'dolphin':'Tursiops_truncatus',
                 'ferret':'Mustela_putorius_furo',
                 'goat':'Capra_hircus',
                 'gorilla':'Gorilla_gorilla_gorilla',
                 'horse':'Equus_caballus',
                 'human':'Homo_sapiens',
                 'mouse':'Mus_musculus',
                 'naked-mole-rat':'Heterocephalus_glaber',
                 'nonhuman-primates':'nonhuman_primates',
                 'pig':'Sus_scrofa',
                 'platypus':'Ornithorhynchus_anatinus',
                 'rabbit':'Oryctolagus_cuniculus',
                 'rat':'Rattus_norvegicus',
                 'salmon':'Salmo_salar',
                 'sheep':'Ovis_aries',
                 'teleostei':'teleostei',
                 'trout':'Oncorhynchus_mykiss',
                 'zebrafish':'Danio_rerio'}

    return species_dict[species]
def _check_format(df, species, chains,verbose, print_invalid_genes = False):
    # map any V genes that can be converted to IMGT.
    # If any additional pairs come up, they can be added to mapper.csv
    mapper = pd.read_csv(DATA_PATH + 'mapper.csv')
    mapper = mapper[mapper['species'] == species]
    ref = pd.read_csv(DATA_PATH + 'imgt_genedb_full.csv')
    ref = ref[ref['species'] == 'Homo sapiens'] #CHANGE THIS PRIOR TO FINAL VERSION
   # ref = ref[ref['species'] == species]

    di = mapper.set_index('v')['imgt'].to_dict()

    df = df.replace({'TRBV': di, 'TRAV': di, 'IGHV': di, 'IGLV': di})
    orig_df = df
    # remove any rows with V genes that are still not in IMGT format
    # these are likely pseudogenes or oddly formatted seqs that haven't been
    # added to the mapper file

    init = df.shape[0]

    if (chains == 'paired'):
        if ('IGHV' and 'IGLV' in df.columns):
            notIn = df['IGHV'][df['IGHV'].isin(ref['imgt']) == False].unique().tolist()
            notIn = notIn + df['IGLV'][df['IGLV'].isin(ref['imgt']) == False].unique().tolist()
            df = df[df['IGHV'].isin(ref['imgt'])]
            df = df[df['IGLV'].isin(ref['imgt'])]
        else:
            notIn = df['TRBV'][df['TRBV'].isin(ref['imgt']) == False].unique().tolist()
            notIn = notIn + df['TRAV'][df['TRAV'].isin(ref['imgt']) == False].unique().tolist()
            df = df[df['TRBV'].isin(ref['imgt'])]
            df = df[df['TRAV'].isin(ref['imgt'])]
    if (chains == 'beta' or chains == 'heavy'):
        if ('IGHV' in df.columns):
            notIn = df['IGHV'][df['IGHV'].isin(ref['imgt']) == False].unique().tolist()
            df = df[df['IGHV'].isin(ref['imgt'])]
        else:
            notIn = df['TRBV'][df['TRBV'].isin(ref['imgt']) == False].astype('string').unique()
            df = df[df['TRBV'].isin(ref['imgt'])]
    if (chains == 'alpha' or chains == 'light'):
        if ('IGLV' in df.columns):
            notIn = df['IGLV'][df['IGLV'].isin(ref['imgt']) == False].unique().tolist()
            df = df[df['IGLV'].isin(ref['imgt'])]
        else:
            notIn = df['TRAV'][df['TRAV'].isin(ref['imgt']) == False].unique().tolist()
            df = df[df['TRAV'].isin(ref['imgt'])]
    print_update('Removed ' + str(init - df.shape[0]) + ' rows not in IMGT format!', verbose)
    if(print_invalid_genes):
        n = '; '.join([str(i) for i in notIn])
        print_update('Unrecognized genes: ' + n, verbose)
    return df, orig_df
def _score_chunks(df, column, metric, mode = 'pairwise', group_column = 'group'):
    # Generate sequence pairs
    seq_list = df[column].values
    seq_pairs = list(itertools.combinations(seq_list, 2)) # no diagonal

    if mode == 'axb': #if mode is pairwise, proceed as normal and skip this step.
       group_list = df[group_column].values
       group_pairs = list(itertools.combinations(group_list, 2))
       ngroups = [ len(set(p)) for p in group_pairs]
       keep = [i for i,j in enumerate(ngroups) if j == 2]
       seq_pairs = [seq_pairs[p] for p in keep]
       id_pairs = list(itertools.combinations(df['Homolig.ID'].values, 2))
       id_pairs = [id_pairs[p] for p in keep]
       id1, id2 = map(list, zip(*id_pairs))
       dim1_labs = sorted(list(set(id1)))
       dim2_labs = sorted(list(set(id2)))
    #Execute scoring
    AA1_sequences, AA2_sequences = map(list,zip(*seq_pairs))
    score_cpp = list(homoligcpp.homolig(DATA_PATH+"align_matrices/"+metric.upper(),  homoligcpp.VectorString(AA1_sequences), homoligcpp.VectorString(AA2_sequences)))
    score_cpp_diag = list(homoligcpp.homolig(DATA_PATH+"align_matrices/"+metric.upper(),  homoligcpp.VectorString(seq_list), homoligcpp.VectorString(seq_list)))
    #Format output matrix
    if mode == 'pairwise':
        # build upper triangle and fill with scores
        tri = np.zeros((len(seq_list),len(seq_list)))
        tri[np.triu_indices(len(seq_list), 1)] = score_cpp
        # add transpose of upper triangle
        tri = tri + tri.T
        # fill diagonal with 1
        np.fill_diagonal(tri,score_cpp_diag)
        score_mat = sparse.csr_matrix(tri)
    if mode == 'axb':
        score_mat = np.reshape(score_cpp, (len(dim1_labs), len(dim2_labs)) )
    return score_mat
def _get_cdrs(fasta_filename):
    CDR1_START = 78
    CDR1_END   = 114
    CDR2_START = 165
    CDR2_END   = 195
    CDR25_START = 240
    CDR25_END = 258
    with open(fasta_filename, 'r') as f:
        sequences = SeqIO.parse(f, 'fasta')
        fasta_with_cdrs = []
        for sequence in sequences:
            # to add full length translated sequence
            #nt_seq = str(sequence.seq).replace('.', '')
            #nt_seq = nt_seq[int(codon_start)-1:]
            #protein_sequence = translate(nt_seq)

            # get data from each record by index
            vals = sequence.description.split('|')
            # do not add any rows with pseudogenes
            #if 'P' in vals[3]:
            #    continue

            # translate nt sequence
            # - only attempt if nt is in-frame. AAG
            nt_sequence = str(sequence.seq)
            cdr1_nt = nt_sequence[CDR1_START:CDR1_END].replace('.', '')
            cdr2_nt = nt_sequence[CDR2_START:CDR2_END].replace('.', '')
            cdr25_nt = nt_sequence[CDR25_START:CDR25_END].replace('.', '')
            if len(cdr1_nt) % 3 == 0:
                cdr1_aa = translate(cdr1_nt)
            else:
                cdr1_aa = '*'
            if len(cdr2_nt) % 3 == 0:
                cdr2_aa = translate(cdr2_nt)
            else:
                 cdr2_aa = '*'
            if len(cdr25_nt) %3 == 0:
                cdr25_aa = translate(cdr25_nt)
            else:
                cdr25_aa = '*'
            fasta_with_cdrs.append(vals[:8] + [cdr1_aa] + [cdr2_aa] + [cdr25_aa])
            vgene_df = pd.DataFrame(fasta_with_cdrs)
            vgene_df.rename(columns={vgene_df.columns[1]: 'V', vgene_df.columns[-3]: 'CDR1', 
                            vgene_df.columns[-2]:'CDR2', vgene_df.columns[-1]:'CDR2.5'}, inplace=True)
            # fill any missing sequences with *
            vgene_df.replace('', '*', inplace=True)
            #vgene_df.to_csv('trv_fasta_cdrs.csv')
    return vgene_df
def _cdr_score(df, c1, c2, c25, metric, mode = 'pairwise', group_column = 'group'):
    # score trv regions
    cdr1_scores = _score_chunks(df, c1, metric, mode, group_column)
    cdr1_scores = cdr1_scores*0.333
    cdr2_scores = _score_chunks(df, c2,  metric, mode, group_column)
    # combine weighted cdr1 and cdr2 score to keep memory usage down
    cdr2_scores = cdr2_scores*0.333 + cdr1_scores
    cdr25_scores = _score_chunks(df, c25, metric, mode, group_column)
    # multiply arrays by 1/3 to scale V gene sequences
    cdr25_scores = cdr2_scores + cdr25_scores*0.333
    return cdr25_scores
def _vgene_msa(df, ref, seq_type, chains, metric,mode = 'pairwise', group_column = 'group'):
    if seq_type == 'tcr' and chains == 'beta':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRBV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5', metric, mode, group_column)
 
    if seq_type == 'tcr' and chains == 'alpha':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRAV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5', metric, mode, group_column)

    if seq_type == 'bcr' and chains == 'heavy':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGHV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5',  metric, mode, group_column)
    if seq_type == 'bcr' and chains == 'light':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGLV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5',metric, mode, group_column)

    if seq_type == 'tcr' and chains == 'paired':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRBV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1b', 'CDR2':'CDR2b', 'CDR2.5':'CDR2.5b'}, inplace=True)
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRAV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1a', 'CDR2':'CDR2a', 'CDR2.5':'CDR2.5a'}, inplace=True)
        alpha_score = _cdr_score(df, 'CDR1a', 'CDR2a', 'CDR2.5a', metric, mode, group_column)
        beta_score = _cdr_score(df, 'CDR1b', 'CDR2b', 'CDR2.5b',  metric, mode, group_column)
        vgene_score = alpha_score + beta_score

    if seq_type == 'bcr' and chains == 'paired':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGHV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1b', 'CDR2':'CDR2b', 'CDR2.5':'CDR2.5b'}, inplace=True)
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGLV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1a', 'CDR2':'CDR2a', 'CDR2.5':'CDR2.5a'}, inplace=True)
        alpha_score = _cdr_score(df, 'CDR1a', 'CDR2a', 'CDR2.5a', metric, mode, group_column)
        beta_score = _cdr_score(df, 'CDR1b', 'CDR2b', 'CDR2.5b', metric, mode, group_column)
        vgene_score = alpha_score + beta_score

    if chains == 'paired':
        alpha_score = _cdr_score(df, 'CDR1a', 'CDR2a', 'CDR2.5a',metric, mode, group_column)
        beta_score = _cdr_score(df, 'CDR1b', 'CDR2b', 'CDR2.5b', metric, mode, group_column)
        vgene_score = alpha_score + beta_score
    #print('vgene array', len(vgene_score)) 
    return df, vgene_score
def _prep_tcr(df, seq_type, chains, metric, species, mode = 'pairwise', group_column = 'group', save_germline=False):
   # get cdr1, cdr2, cdr2.5 sequences from trv names
    tra = DATA_PATH + 'fastas/' + species + '/TR/TRAV.fasta'
    trb = DATA_PATH + 'fastas/' + species + '/TR/TRBV.fasta'
    tra_df = _get_cdrs(tra)
    trb_df = _get_cdrs(trb)
    trv_df = pd.concat([tra_df, trb_df])
    # get score dictionaries for trv cdrs
    df,  germline_scores  = _vgene_msa(df, trv_df, seq_type, chains,  metric, mode, group_column)
    return df, germline_scores
def _prep_bcr(df, seq_type, chains, metric, species, mode = 'pairwise', group_column = 'group'):
    # get cdr1, cdr2, cdr2.5 sequences from igv names
    igh = DATA_PATH + 'fastas/' + species + '/IG/IGHV.fasta'
    igl = DATA_PATH + 'fastas/' + species + '/IG/IGLV.fasta'
    igk = DATA_PATH + 'fastas/' + species + '/IG/IGKV.fasta'
    igh_df = _get_cdrs(igh)
    igl_df = _get_cdrs(igl)
    igk_df = _get_cdrs(igk)
    igv_df = pd.concat([igh_df, igl_df, igk_df])
    # get score dictionaries for _trv cdrs
    df, germline_scores = _vgene_msa(df, igv_df, seq_type, chains, metric, mode, group_column)
    return df, germline_scores
def _consolidate_input(df, colnames,verbose):
	#Collapse dataframe to only contain unique sequences (remove duplicates) based on fields provided in colnames vector
	#Label every sequence with a unique identifier.
	#Print result (# unique seqs removed) and return consolidated input.
    pasted = []

    for r in range(df.shape[0]):
        pasted.append('')

    for l in range(len(colnames)):
        pasted = pasted + df[colnames[l]].astype(str)

    pasted = pasted.tolist()
    #If not present already, add a 'Count' column in data.
    #Column must be prenamed to exactly 'Counts' to be recognized.
    #Using pre-defined counts assumes sequences are already unique, or else that the provided
    #counts are the sum of all non-unique sequences already.
    if 'Count' not in df.columns:
        counts = [pasted.count(x) for x in pasted]
        df['Count'] = counts

    unique_indices = sorted([pasted.index(x) for x in set(pasted)])
    collapsed = df[colnames].iloc[unique_indices,:]
    collapsed['Count'] = df['Count'].iloc[unique_indices]
    # Make Homolig ID column, which is guaranteed to be a string.
    # prevents ImplicitModificationWarning converting to AnnData later.
    lst = range(collapsed.shape[0])
    lst = [format(x, 'd') for x in lst]
    collapsed['Homolig.ID'] = ['H' + s for s in lst]

    n_removed = df.shape[0] - collapsed.shape[0]
    print_update('Removed '+ str(n_removed) + ' non-unique sequences!', verbose)
    collapsed = df #totally nullifying this function. Can re-institute if needed.
    print_update(str(collapsed.shape[0]) + ' sequences remaining.', verbose)

    return collapsed
def _prep_input(input_file, seq_type, chains, metric, species, mode, verbose, group_column = 'group', print_invalid_genes = False):
    # read in input file
    df = pd.read_csv(input_file, dtype='category')
    df, orig_df = _check_format(df, species, chains, verbose, print_invalid_genes=print_invalid_genes)

    #Collapse input to only include relevant columns and unique sequences.
    if seq_type == 'bcr':
        if chains == 'heavy':
            df = _consolidate_input(df, ['IGHV', 'CDR3.heavy.aa'],verbose)
        if chains == 'light':
            df = _consolidate_input(df, ['IGLV', 'CDR3.light.aa'],verbose)
        if chains == 'paired':
            df = _consolidate_input(df, ['IGHV', 'CDR3.heavy.aa', 'IGLV', 'CDR3.light.aa'],verbose)
    if seq_type == 'tcr':
        if chains == 'alpha':
            df = _consolidate_input(df, ['TRAV', 'CDR3.alpha.aa'],verbose)
        if chains == 'beta':
            df = _consolidate_input(df, ['TRBV', 'CDR3.beta.aa'],verbose)
        if chains == 'paired':
            df = _consolidate_input(df, ['TRAV', 'CDR3.alpha.aa', 'TRBV', 'CDR3.beta.aa'],verbose)
    if seq_type == 'seq':
        df = _consolidate_input(df, ['seq'],verbose)

    df['Homolig.ID'] = df.index
    return df, orig_df
def _make_output(data, df, mode = 'pairwise', group_column = 'group'):
    #df.Count = [int(i) for i in df.Count]
    df.Count = [int(i) if not pd.isna(i) else -1 for i in df.Count]
    if mode == 'pairwise':
        adata = ad.AnnData(data, obs = df, dtype = 'float32')
    if mode == 'axb':
        groups = sorted(list(set(df[group_column])))
        obs_df = df[df[group_column] == groups[0]]
        var_df =  df[df[group_column] == groups[1]]
        adata = ad.AnnData(data, obs = obs_df, var = var_df, dtype = 'float32')
    return adata
def _default_output_filename(input_file, seq_type, chains, mode = 'pairwise'):
    #Return the appropriate default output file name for anndata object.
    #Doing within this function such that outputs can be handled within master homolig function rather than
    #chain-specific subfunctions.

    if mode == 'axb':
        modesuffix='_axb'
    else:
        modesuffix=''

    if seq_type == 'tcr' and chains == 'beta':
        filename = os.path.splitext(input_file)[0] + modesuffix + "_beta_tcr.h5ad"
    elif seq_type == 'tcr' and chains == 'alpha':
        filename = os.path.splitext(input_file)[0] + modesuffix +  "_alpha_tcr.h5ad"
    elif seq_type == 'tcr' and chains == 'paired':
        filename = os.path.splitext(input_file)[0] +  modesuffix + "_paired_tcr.h5ad"
    elif seq_type == 'bcr' and chains == 'heavy':
        filename = os.path.splitext(input_file)[0] +  modesuffix +  "_heavy_bcr.h5ad"
    elif seq_type == 'bcr' and chains == 'light':
        filename = os.path.splitext(input_file)[0] +  modesuffix + "_light_bcr.h5ad"
    elif seq_type == 'bcr' and chains == 'paired':
        filename = os.path.splitext(input_file)[0] +  modesuffix + "_paired_bcr.h5ad"
    elif seq_type == 'seq':
        filename = os.path.splitext(input_file)[0] +  modesuffix +  "_seq.h5ad"

    return filename
def print_update(message, verbose):
    #Get current time. Print to console along w message, if verbose==True
    if verbose:
        now = strftime("%Y-%m-%d %H:%M:%S", localtime())
        print('[' + now + ']', message)

def homolig_format_checker(input_file, seq_type, chains=None, metric='aadist', species='human', mode='pairwise', input2 = None,
            output_file = None, verbose = True, save_germline = False, save_reformatted_input = False):
    # Print start-up message
    now = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('[' + now + ']', 'Homolig version 0.1 Format Checker')  # Find a way to add __version__ attribute to package at later date.
    print('                         input_file: ', input_file)
    print('                           seq_type: ', seq_type)
    print('                             chains: ', chains)
    print('                             metric: ', metric)
    print('                            species: ', species)
    print('                               mode: ', mode)
    print('                            verbose: ', verbose)
    if(mode == 'axb'):
        print('                             input2:', input2)
    if output_file is None:
        output_file = _default_output_filename(input_file, seq_type, chains, mode)
    print('                             output: ', output_file)
    if save_germline:
        print('                      save_germline:', save_germline)

    if save_reformatted_input:
        print('             save_reformatted_input:', save_reformatted_input)


    valid_seq_types = ['tcr', 'bcr', 'seq']
    valid_chains = ['alpha', 'beta', 'paired', 'light', 'heavy', None]
    valid_metrics = ['aadist', 'aadist_euc', 'blosum62']
    valid_species = ['human', 'mouse', 'mas-night-monkey', 'rhesus-monkey',
                     'alpaca', 'bovine', 'camel', 'catfish', 'chicken', 'chondrichthyes'
                                                                        'cod', 'crab-eating-macaque', 'dog', 'dolphin',
                     'ferret',
                     'goat', 'gorilla', 'horse', 'naked-mole-rat',
                     'nonhuman-primates', 'pig', 'platypus', 'rabbit', 'rat',
                     'salmon', 'sheep', 'teleostei', 'trout', 'zebrafish']
    valid_modes = ['pairwise', 'axb']
    if seq_type not in valid_seq_types:
        raise ValueError('Not valid sequence type')
    if seq_type == 'seq':
        chains = None
    if chains not in valid_chains:
        raise ValueError('Not valid chain type')
    # if metric not in valid_metrics:
    #    raise ValueError('Not valid sequence type')
    if species not in valid_species:
        raise ValueError('Not valid species type')
    if mode not in valid_modes:
        raise ValueError('Not valid mode')

    #Convert species to latin before continuing: database files are annotated this way.
    #This is needlessly complex and should be fixed eventually.
    species = _species_lookup(species)
    #Consolidate input prior to continuing.
    if mode == 'pairwise':
        group_column = 'group'#placeholder dummy variable to pass into sub-functions.
        df, orig_df  = _prep_input(input_file, seq_type, chains, metric, species, mode,verbose, print_invalid_genes=True)
    #If axb mode, load second input file and concatenate into one input.
    if mode == 'axb':
        print_update('Processing input 1:',verbose)
        df, orig_df = _prep_input(input_file, seq_type, chains, metric, species, mode,verbose, print_invalid_genes=True)
        print_update('Processing input 2:',verbose)
        df2, orig_df2 = _prep_input(input2, seq_type, chains, metric, species, mode,verbose, print_invalid_genes=True)
        df['group'] = 'a'
        df2['group'] = 'b'
        df = pd.concat([df,df2])
        orig_df = pd.concat[orig_df, orig_df2]
        lst = range(df.shape[0])
        lst = [format(x, 'd') for x in lst]
        df['Homolig.ID'] = ['H' + s for s in lst]

    if save_reformatted_input:
        reformatted_input_file = output_file[:-5] + '_reformatted-input.csv'
        orig_df.to_csv(reformatted_input_file, index=False)
        abridged_input_file = output_file[:-5] + '_abridged-input.csv'
        df.to_csv(abridged_input_file, index = False)

    now = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('[' + now + ']','Homolig Format Checker completed.')
    return 0
def homolig(input_file, seq_type, chains=None, metric='aadist', species='human', mode='pairwise', input2 = None,
            output_file=None, verbose=False, save_germline=False, save_reformatted_input=False):
    # Print start-up message
    now = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('[' + now + ']', 'Homolig version 0.1')  # Find a way to add __version__ attribute to package at later date.
    print('                         input_file: ', input_file)
    print('                           seq_type: ', seq_type)
    print('                             chains: ', chains)
    print('                             metric: ', metric)
    print('                            species: ', species)
    print('                               mode: ', mode)
    print('                            verbose: ', verbose)
    if(mode == 'axb'):
        print('                             input2:', input2)
    if output_file is None:
        output_file = _default_output_filename(input_file, seq_type, chains, mode)
    print('                             output: ', output_file)
    if save_germline:
        print('                      save_germline:', save_germline)
    if save_reformatted_input:
        print('             save_reformatted_input:', save_reformatted_input)

    valid_seq_types = ['tcr', 'bcr', 'seq']
    valid_chains = ['alpha', 'beta', 'paired', 'light', 'heavy', None]
    valid_metrics = ['aadist', 'aadist_euc', 'blosum62']
    valid_species = ['human', 'mouse', 'mas-night-monkey', 'rhesus-monkey',
                     'alpaca', 'bovine', 'camel', 'catfish', 'chicken', 'chondrichthyes'
                                                                        'cod', 'crab-eating-macaque', 'dog', 'dolphin',
                     'ferret',
                     'goat', 'gorilla', 'horse', 'naked-mole-rat',
                     'nonhuman-primates', 'pig', 'platypus', 'rabbit', 'rat',
                     'salmon', 'sheep', 'teleostei', 'trout', 'zebrafish']
    valid_modes = ['pairwise', 'axb']
    if seq_type not in valid_seq_types:
        raise ValueError('Not valid sequence type')
    if seq_type == 'seq':
        chains = None
    if chains not in valid_chains:
        raise ValueError('Not valid chain type')
    # if metric not in valid_metrics:
    #    raise ValueError('Not valid sequence type')
    if species not in valid_species:
        raise ValueError('Not valid species type')
    if mode not in valid_modes:
        raise ValueError('Not valid mode')

    #Convert species to latin before continuing: database files are annotated this way.
    #This is needlessly complex and should be fixed eventually.
    species = _species_lookup(species)
    #Consolidate input prior to continuing.
    if mode == 'pairwise':
        group_column = 'group'#placeholder dummy variable to pass into sub-functions.
        df, orig_df  = _prep_input(input_file, seq_type, chains, metric, species, mode,verbose)
    #If axb mode, load second input file and concatenate into one input.
    if mode == 'axb':
        print_update('Processing input 1:',verbose)
        df, orig_df = _prep_input(input_file, seq_type, chains, metric, species, mode,verbose)
        print_update('Processing input 2:',verbose)
        df2, orig_df2 = _prep_input(input2, seq_type, chains, metric, species, mode,verbose)
        df['group'] = 'a'
        df2['group'] = 'b'
        df = pd.concat([df,df2])
        lst = range(df.shape[0])
        lst = [format(x, 'd') for x in lst]
        df['Homolig.ID'] = ['H' + s for s in lst]

    if seq_type == 'tcr':
        # preprocessing
        print_update('scoring germline V regions', verbose)
        df, germline_scores = _prep_tcr(df, seq_type, chains, metric, species, mode)

        if chains == 'alpha' or chains == 'paired':
            print_update('scoring CDR3a',verbose)
            cdr3a_scores = _score_chunks(df, 'CDR3.alpha.aa', metric, mode)
        if chains == 'beta' or chains == 'paired':
            print_update('scoring CDR3b', verbose)
            cdr3b_scores = _score_chunks(df, 'CDR3.beta.aa', metric, mode)

        print_update('merging dictionaries', verbose)
        if chains == 'beta':
            all_scores = germline_scores + cdr3b_scores
        elif chains == 'alpha':
            all_scores = germline_scores + cdr3a_scores
        elif chains == 'paired':
            all_scores = germline_scores + cdr3a_scores + cdr3b_scores

        adata = _make_output(all_scores, df, mode=mode)

    elif seq_type == 'seq':
        print_update('scoring sequences', verbose)
        seq_scores = _score_chunks(df, 'seq', metric, mode, group_column)
        adata = _make_output(seq_scores, df, mode, group_column)

    elif seq_type == 'bcr':
        # preprocessing
        print_update('scoring germline V regions', verbose)
        df, germline_scores = _prep_bcr(df, seq_type, chains, metric, species, mode)

        if chains == 'light' or chains == 'paired':
            print_update('scoring CDR3l',verbose)
            cdr3l_scores = _score_chunks(df, 'CDR3.light.aa', metric, mode)
        if chains == 'heavy' or chains == 'paired':
            print_update('scoring CDR3b', verbose)
            cdr3h_scores = _score_chunks(df, 'CDR3.heavy.aa', metric, mode)

        print_update('merging dictionaries', verbose)
        if chains == 'heavy':
            all_scores = germline_scores + cdr3h_scores
        elif chains == 'light':
            all_scores = germline_scores + cdr3l_scores
        elif chains == 'paired':
            all_scores = germline_scores + cdr3h_scores + cdr3l_scores

        adata = _make_output(all_scores, df, mode=mode)

    else:
        raise ValueError('Not valid sequence type')

    adata.write(output_file)

    if save_germline:
        germ_file = output_file[:-5] + '_germline.h5ad'
        adata_germ = _make_output(germline_scores, df, mode=mode)
        adata_germ.write(germ_file)

        for cdr3dat in ['cdr3b_scores', 'cdr3a_scores', 'cdr3h_scores', 'cdr3l_scores']:
            if cdr3dat in locals():
                adata_cdr3 = _make_output(locals()[cdr3dat], df, mode = mode)
                cdr3_file = output_file[:-5] + '_' + cdr3dat[:-7] + '.h5ad'
                adata_cdr3.write(cdr3_file)
    if save_reformatted_input:
        reformatted_input_file = output_file[:-5] + '_reformatted-input.csv'
        orig_df.to_csv(reformatted_input_file, index=False)

    now = strftime("%Y-%m-%d %H:%M:%S", localtime())
    print('[' + now + ']','Homolig completed.')
    return 0
