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

def _check_format(df, species, chains):
    # map any V genes that can be converted to IMGT.
    # If any additional pairs come up, they can be added to mapper.csv
    d = pd.read_csv(DATA_PATH + '/mapper.csv')
    d = d[d['species'] == species]
    di = d.set_index('v')['imgt'].to_dict()
    df = df.replace({'TRBV': di, 'TRAV': di, 'IGHV': di, 'IGLV': di})
    # remove any rows with V genes that are still not in IMGT format
    # these are likely pseudogenes or oddly formatted seqs that haven't been
    # added to the mapper file
    ref = pd.read_csv(DATA_PATH + '/imgt_genedb_full.csv')
    ref = ref[ref['species'] == species]
    init = df.shape[0]
    if (chains == 'paired'):
        if ('IGHV' and 'IGLV' in df.columns):
            df = df[~df['IGHV'].isin(ref['imgt'])]
            df = df[~df['IGLV'].isin(ref['imgt'])]
        else:
            df = df[~df['TRBV'].isin(ref['imgt'])]
            df = df[~df['TRAV'].isin(ref['imgt'])]
    if (chains == 'beta' or chains == 'heavy'):
        if ('IGHV' in df.columns):
            df = df[~df['IGHV'].isin(ref['imgt'])]
        else:
            df = df[~df['TRBV'].isin(ref['imgt'])]
    if (chains == 'alpha' or chains == 'light'):
        if ('IGLV' in df.columns):
            df = df[~df['IGLV'].isin(ref['imgt'])]
        else:
            df = df[~df['TRAV'].isin(ref['imgt'])]
    print('Removed ', init - df.shape[0], 'rows not in IMGT format!')

    # Make Homolig ID column, which is guaranteed to be a string.
    # prevents ImplicitModificationWarning converting to AnnData later. AAG
    lst = range(df.shape[0])
    lst = [format(x,'d') for x in lst]
    df['Homolig.ID'] = ['H'+ s for s in lst ]
    df = df.set_index('Homolig.ID')
    return df

def _get_matrix(metric):
    matpath = DATA_PATH + '/align_matrices/' + metric.upper()
    with open(matpath) as handle:
        #print(matrix.alphabet)
        #print(matrix['A','D'])
        return substitution_matrices.read(handle)


# pairwise scoring to dictionary 
def _score_pairwise(aln_list, key_list, score_matrix):
    keys = []
    values = []
    # initialize score 
    score = 0
    for index, aln in enumerate(list(aln_list)):
    #seq_pairs = list(itertools.combinations(aln_list, 2))
    #for seq_pair in seq_pairs:
        #print(aln)
        seq1 = aln[0]
        seq2 = aln[1]
        # length of amino acid characters
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        # get list of sliding windows based on min length
        s1 = ["".join(w) for w in mit.windowed(seq1, min_len)]
        s2 = ["".join(w) for w in mit.windowed(seq2, min_len)]
        # get tuples of sequence comparison pairs
        pairs = [list(tup) for tup in itertools.product(s1, s2)]
        # convert list of strings to lists of aa characters
        # get tuple of aligned aas in each sequence to compare)
        max_score = 0
        for pair in pairs:
            #print(pair)
            if pair[0] == pair[1]:
                max_score = 1.0
                break
            aa_pairs = list(zip(pair[0], pair[1]))
            result = sum(score_matrix[aa_pair] for aa_pair in aa_pairs) / max_len
            if result > max_score:
                max_score = result
        keys += [tuple(sorted(aln))]
        #keys += [tuple(sorted(list(aln_list)[index]))]
        values += [max_score]
    # returns tuple of ((sequnce tuple): score)
    #pprint(list(zip(keys, values)))
    #return zip(keys,values)
    return values

#def chunks(iterable, chunk_size=200):
#    iterator = iter(iterable)
#    for first in iterator:
#        yield itertools.chain([first], itertools.islice(iterator, chunk_size - 1))

def _score_chunks(df, column, score_matrix, metric, chunk_size=200):
    # get list of unique sequences 
    #df = df.drop_duplicates(subset=[column])
    #seq_list = df[column].tolist()
    seq_list = df[column].values
    seq_pairs = list(itertools.combinations(seq_list, 2)) # no diagonal
    AA1_sequences, AA2_sequences = map(list,zip(*seq_pairs))
    # seq_chunks = [seq_pairs[i:i + chunk_size] for i in range(0, len(seq_pairs), chunk_size)]
    # zip_args = list(zip(seq_chunks, seq_chunks))
    # chunk_args = [x + (score_matrix,) for x in zip_args]
    # with Pool(cpu_count()) as pool:
    #     chunk_scores = pool.starmap(_score_pairwise, chunk_args)

    # score = list(itertools.chain(*chunk_scores))
    score_cpp = list(homoligcpp.homolig(DATA_PATH+"align_matrices/"+metric.upper(),  homoligcpp.VectorString(AA1_sequences), homoligcpp.VectorString(AA2_sequences)))
    # build upper triangle and fill with scores 
    tri = np.zeros((len(seq_list),len(seq_list)))
    tri[np.triu_indices(len(seq_list), 1)] = score
    # add transpose of upper triangle
    tri = tri + tri.T
    # fill diagonal with 1 
    np.fill_diagonal(tri,1)
    score_mat = sparse.csr_matrix(tri)
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

def _cdr_score(df, c1, c2, c25, score_matrix, metric):
    # score trv regions
    cdr1_scores = _score_chunks(df, c1, score_matrix, metric)
    cdr1_scores = cdr1_scores*0.333
    cdr2_scores = _score_chunks(df, c2, score_matrix, metric)
    # combine weighted cdr1 and cdr2 score to keep memory usage down
    cdr2_scores = cdr2_scores*0.333 + cdr1_scores
    cdr25_scores = _score_chunks(df, c25, score_matrix, metric)
    # multiply arrays by 1/3 to scale V gene sequences
    cdr25_scores = cdr2_scores + cdr25_scores*0.333
    return cdr25_scores

def _vgene_msa(df, ref, seq_type, chains, metric, score_matrix):
    if seq_type == 'tcr' and chains == 'beta':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRBV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5', score_matrix, metric)
 
    if seq_type == 'tcr' and chains == 'alpha':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRAV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5', score_matrix, metric)

    if seq_type == 'bcr' and chains == 'heavy':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGHV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5', score_matrix, metric)
    if seq_type == 'bcr' and chains == 'light':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGLV'], right_on= ['V'], how = 'left')
        vgene_score = _cdr_score(df, 'CDR1', 'CDR2', 'CDR2.5', score_matrix, metric)

    if seq_type == 'tcr' and chains == 'paired':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRBV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1b', 'CDR2':'CDR2b', 'CDR2.5':'CDR2.5b'}, inplace=True)
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['TRAV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1a', 'CDR2':'CDR2a', 'CDR2.5':'CDR2.5a'}, inplace=True)
        alpha_score = _cdr_score(df, 'CDR1a', 'CDR2a', 'CDR2.5a', score_matrix, metric)
        beta_score = _cdr_score(df, 'CDR1b', 'CDR2b', 'CDR2.5b', score_matrix, metric)
        vgene_score = alpha_score + beta_score

    if seq_type == 'bcr' and chains == 'paired':
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGHV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1b', 'CDR2':'CDR2b', 'CDR2.5':'CDR2.5b'}, inplace=True)
        df = pd.merge(df, ref[['V','CDR1','CDR2','CDR2.5']], left_on=  ['IGLV'], right_on= ['V'], how = 'left')
        df.rename(columns={'CDR1':'CDR1a', 'CDR2':'CDR2a', 'CDR2.5':'CDR2.5a'}, inplace=True)
        alpha_score = _cdr_score(df, 'CDR1a', 'CDR2a', 'CDR2.5a', score_matrix, metric)
        beta_score = _cdr_score(df, 'CDR1b', 'CDR2b', 'CDR2.5b', score_matrix, metric)
        vgene_score = alpha_score + beta_score

    #print('vgene array', len(vgene_score)) 
    return df, vgene_score

def homolig(input_file, seq_type, chains=None, metric='aadist', species='human'):
    # validate_input(**locals())
    valid_seq_types = ['tcr', 'bcr', 'seq']
    valid_chains = ['alpha', 'beta', 'paired', 'light', 'heavy', None]
    valid_metrics = ['aadist', 'aadist_euc', 'blosum62']
    valid_species = ['human', 'mouse', 'mas-night-monkey', 'rhesus-monkey',
                     'alpaca', 'bovine', 'camel', 'catfish', 'chicken', 'chondrichthyes'
                     'cod', 'crab-eating-macaque', 'dog', 'dolphin', 'ferret',
                     'goat', 'gorilla', 'horse', 'naked-mole-rat',
                     'nonhuman-primates', 'pig', 'platypus', 'rabbit', 'rat',
                     'salmon', 'sheep', 'teleostei', 'trout', 'zebrafish']   
    if seq_type not in valid_seq_types:
        raise ValueError('Not valid sequence type')
    if chains not in valid_chains:
        raise ValueError('Not valid chain type')
    if metric not in valid_metrics:
        raise ValueError('Not valid sequence type')
    if species not in valid_species:
        raise ValueError('Not valid species type')

    if seq_type == 'tcr':
        # preprocessing
        df, score_matrix, trv_aligned, germline_scores = _prep_tcr(input_file, seq_type, chains,  metric, species)
        if chains == 'paired':
            adata = _tcr_paired(input_file, metric, df, score_matrix, germline_scores)
        elif chains == 'alpha':
            adata = _tcr_alpha(df,  score_matrix, germline_scores, input_file, metric)
        elif chains == 'beta':
            #if any(~df['TRBV'].isin(trv_aligned['V'])):
            #    raise ValueError('V genes are not in valid IMGT format')
            #else:
            adata = _tcr_beta(df, score_matrix, germline_scores, input_file, metric)
        else:
            raise ValueError('Not valid sequence type')

    elif seq_type == 'seq':
        adata = _seq(input_file, metric)
    
    elif seq_type == 'bcr':
        # preprocessing
        df, igv_aligned, score_matrix, germline_scores  = _prep_bcr(input_file, seq_type, chains,  metric, species)
        #print(igv_aligned)
        if chains == 'paired':
            # check v gene format
            adata = _bcr_paired(input_file, metric, df, score_matrix, germline_scores)
        elif chains == 'light':
            adata = _bcr_light(df, score_matrix, germline_scores, input_file, metric)
        elif chains == 'heavy':
            adata = _bcr_heavy(df, score_matrix, germline_scores, input_file, metric)

    return adata

def _prep_tcr(input_file, seq_type, chains, metric, species): 
    # read in input file 
    df = pd.read_csv(input_file, dtype='category')
    df = _check_format(df, species, chains)
    species = _species_lookup(species)
    #print(df.memory_usage(deep=True))
    #print(df.dtypes)
    # read in score matrix
    score_matrix = _get_matrix(metric)
    # get cdr1, cdr2, cdr2.5 sequences from trv names
    tra = DATA_PATH + '/fastas/' + species + '/TR/TRAV.fasta'
    trb = DATA_PATH + '/fastas/' + species + '/TR/TRBV.fasta'
    tra_df = _get_cdrs(tra)
    trb_df = _get_cdrs(trb)
    trv_df = pd.concat([tra_df, trb_df])
    #trv_df.to_csv('trv_test.csv')
    # get score dictionaries for trv cdrs
    trv_aligned, germline_scores = _vgene_msa(df, trv_df, seq_type, chains,  metric, score_matrix)
    return df, score_matrix, trv_aligned, germline_scores 

def _tcr_paired(input_file, metric, df, score_matrix, germline_scores):
    print('scoring CDR3a')
    cdr3a_scores = _score_chunks(df, 'CDR3.alpha.aa', score_matrix, metric)
    print('scoring CDR3b')
    cdr3b_scores = _score_chunks(df, 'CDR3.beta.aa', score_matrix, metric)
    print('merging dictionaries') 
    # combining key and value lists
    all_scores = germline_scores + cdr3a_scores + cdr3b_scores

    filename = os.path.splitext(input_file)[0] + "_paired_tcr.h5ad"
    adata = ad.AnnData(all_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata

def _tcr_alpha(df, score_matrix, germline_scores, input_file,  metric):
    print('scoring CDR3a')
    cdr3a_scores = _score_chunks(df, 'CDR3.alpha.aa', score_matrix, metric)
    print('merging dictionaries')
    # combining key and value lists
    all_scores = germline_scores + cdr3a_scores
    filename = os.path.splitext(input_file)[0] + "_alpha_tcr.h5ad"
    adata = ad.AnnData(all_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata


def _tcr_beta(df, score_matrix, germline_scores, input_file, metric):
    print('scoring CDR3b')
    cdr3b_scores = _score_chunks(df, 'CDR3.beta.aa', score_matrix, metric)
    print('merging dictionaries')
    # combining key and value lists
    all_scores = germline_scores + cdr3b_scores
    #all_scores = all_scores.reshape(len(df), len(df))
    # arrays of scores
    filename = os.path.splitext(input_file)[0] + "_beta_tcr.h5ad"
    adata = ad.AnnData(all_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata

def _seq(input_file, metric): 
    # read in input file 
    df = pd.read_csv(input_file, dtype='category')
    df = df.astype(str)
    # Make Homolig ID column, which is guaranteed to be a string.
    # prevents ImplicitModificationWarning converting to AnnData later. AAG
    lst = range(df.shape[0])
    lst = [format(x, 'd') for x in lst]
    df['Homolig.ID'] = ['H' + s for s in lst]
    df = df.set_index('Homolig.ID')
    # read in score matrix
    score_matrix = _get_matrix(metric)
    seq_scores = _score_chunks(df, 'sequence', score_matrix, metric)
    filename = os.path.splitext(input_file)[0] + "_sequence.h5ad"
    adata = ad.AnnData(seq_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata

def _prep_bcr(input_file, seq_type, chains, metric, species): 
    # read in input file 
    df = pd.read_csv(input_file, dtype='category')
    df = _check_format(df, species, chains)
    species = _species_lookup(species)
    # read in score matrix
    score_matrix = _get_matrix(metric)
    # get cdr1, cdr2, cdr2.5 sequences from igv names
    igh = DATA_PATH + 'fastas/' + species + '/IG/IGHV.fasta'
    igl = DATA_PATH + 'fastas/' + species + '/IG/IGLV.fasta'
    igk = DATA_PATH + 'fastas/' + species + '/IG/IGKV.fasta'
    igh_df = _get_cdrs(igh)
    igl_df = _get_cdrs(igl)
    igk_df = _get_cdrs(igk)
    
    igv_df = pd.concat([igh_df, igl_df, igk_df])
    # get score dictionaries for _trv cdrs
    igv_aligned, germline_dict = _vgene_msa(df, igv_df, seq_type, chains, metric, score_matrix)
    
    return df, igv_aligned, score_matrix, germline_dict 

def _bcr_paired(input_file, metric, df, score_matrix, germline_scores):
    print('scoring CDR3l')
    cdr3l_scores = _score_chunks(df, 'CDR3.light.aa', score_matrix, metric)
    print('scoring CDR3h')
    cdr3h_scores = _score_chunks(df, 'CDR3.heavy.aa', score_matrix, metric)
    print('merging dictionaries')
    all_scores = germline_scores + cdr3l_scores + cdr3h_scores
    filename = os.path.splitext(input_file)[0] + "_paired_bcr.h5ad"
    adata = ad.AnnData(all_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata

def _bcr_light(df, score_matrix, germline_scores, input_file,  metric):
    print('scoring CDR3l')
    cdr3l_scores = _score_chunks(df, 'CDR3.light.aa', score_matrix, metric)
    print('merging dictionaries')
    # arrays of scores
    all_scores = germline_scores + cdr3l_scores
    filename = os.path.splitext(input_file)[0] + "_light_bcr.h5ad"
    adata = ad.AnnData(all_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata

def _bcr_heavy(df, score_matrix, germline_scores, input_file, metric):
    print('scoring CDR3h')
    cdr3h_scores = _score_chunks(df, 'CDR3.heavy.aa', score_matrix, metric)
    print('merging dictionaries')
    all_scores = germline_scores + cdr3h_scores
    # arrays of scores
    filename = os.path.splitext(input_file)[0] + "_heavy_bcr.h5ad"
    adata = ad.AnnData(all_scores, obs=df, dtype='float32')
    adata.write(filename)
    return adata
