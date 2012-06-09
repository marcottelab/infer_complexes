import examples as ex
import complex as co
import seqs
import orth
import utils as ut
import elution as el
from Struct import Struct
import pairarray as pa
import numpy as np
import fnet
import score
from numpy import zeros,ndarray

def learning_examples(species, seqdb, elut_fs, scores=['poisson','wcc','apex'],
                      fnet_file=None, splits=[0,.33,.66,1], neg_ratios=[10,40],
                      ind_cycle=[0,-1], cutoff=0.25, base_pairdict=None,
                      pos_splits=None):
    """
    Species: 'Hs', 'Ce', ...
    seqdb: 'uni', 'ensp'
    Use fnet_file = -1 to skip functional network. Use default for it to be
      figured out from species and seqdb as in ut.config.
    Provide base_pairdict as [extr,exte] to use saved examples.
    """
    elut_prots = load_prot_set(elut_fs) # for test set negatives
    if base_pairdict is None:
        ppi_cxs,clean_cxs = load_training_complexes(species, seqdb)
        train,test = ex.base_examples(ppi_cxs, clean_cxs, elut_prots, splits,
            nratio_train=neg_ratios[0], nratio_test=neg_ratios[1],
            pos_splits=pos_splits, ind_cycle=ind_cycle)
    else:
        train,test = base_pairdict
    ntest_pos = len([v for k,v in test.d.items() if v[0]==1])
    print 'total test positives:', ntest_pos
    train_test = [score_and_filter(pdict, scores, elut_fs, cutoff, species,
                      seqdb, fnet_file) for pdict in [train,test]]
    print 'done.', stats(train_test)
    return train_test, ntest_pos

def score_and_filter(pdict, scores, elut_fs, cutoff, species, seqdb, fnet_file):
    arr = new_score_array(pdict, scores, elut_fs, fnet_file)
    score.scores_array(arr, elut_fs, scores, cutoff)
    print 'filtering.'
    if cutoff != -1:
        columns = range(3,len(arr[0]))
        # note that this double-filters often. only an efficiency issue.
        arr = filter_arr(arr, columns, cutoff)
    if fnet_file!=-1:
        fnet.score_arr(arr, species, seqdb, fnet_file)
    return arr 

def new_score_array(pd, scores, fs, fnet_file):
    """
    arr[0] returns the first row of data.
    row[name_score(f,s)] returns the proper score
    row['id1'],'id2','hit' index to the first three columns
    """
    data_names = [score.name_score(f,s) for f in fs for s in scores]
    if fnet_file != -1:
        data_names += fnet.fnet_names(fnet_file)
    arr = zeros(shape = len(pd.d),
                  dtype = ','.join(['a20']*2 + ['i1'] + ['f4']*len(data_names)))
    names = ['id1','id2','hit'] + data_names
    arr.dtype.names = tuple(names)
    for i,((p1,p2),lhit) in enumerate(pd.d.items()):
        row = arr[i]
        row['id1'],row['id2'],row['hit'] = p1,p2,lhit
    return arr

def predict_all(species, seqdb, elut_fs, scores=['poisson','wcc','apex'],
                fnet_file=None, score_cutoff=0.25):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    pairs = el.all_filtered_pairs(elut_fs, scores, score_cutoff)
    print len(pairs), 'total interactions passing cutoff'
    exs = [[p1, p2, '?'] for p1,p2 in pairs]
    ex_struct = Struct(examples=exs,names=['id1','id2','hit'])
    el.score_multi_exs([ex_struct], elut_fs, scores, score_cutoff)
    if fnet_file!=-1:
        fnet.score_examples(ex_struct, species, seqdb, fnet_file)
    return ex_struct

def load_prot_set(elut_fs):
    return reduce(set.union, (set(el.load_elution(f).prots) for f in elut_fs))

def stats(train_test):
    nums =  [len([1 for row in arr if row['hit']==tf]) for arr in train_test for tf in [1,0]]
    return 'train %sP/%sN; test %sP/%sN' % tuple(nums)
    
def load_training_complexes(species, seqdb):
    ppi = co.load_complexes_singleline(ut.proj_path('ppi_cxs'))
    clean = co.load_complexes_multiline(ut.proj_path('clean_cxs'))
    ppi,clean = [convert_complexes(cxs, species, seqdb) for cxs in [ppi,clean]]
    return ppi,clean

def convert_complexes(cxs, species, seqdb):
    #assert species=='Hs' or seqdb=='ensp', "translation not supported"
    if seqdb=='ensp':
        dconv = ut.load_dict_sets(ut.proj_path('convert','Hs_uni2Hs_ensp.tab'))
        ps = seqs.load_prots_from_fasta(ut.proj_path('ensp_fasta','Hs_longest'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    if species!='Hs':
        from_sp = 'Hs_'+seqdb
        to_sp = species+'_'+seqdb
        dconv = orth.odict(from_sp, to_sp)
        ps = seqs.load_prots_from_fasta(ut.proj_path('ensp_fasta',
            species+'_longest'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    return cxs

def filter_arr(arr, columns, cutoff):
    """
    Uses max: return a new array of all rows from input arr with at least one
    value in columns exceeding cutoff.
    """
    newarr = arr[[i for i in range(len(arr)) if max([x for j,x in
                enumerate(arr[i]) if j in columns]) > cutoff]]
    del arr
    return newarr
    
def preds_thresh(tested, thresh):
    limit = None
    for i,t in enumerate(tested):
        if t[2]<thresh:
            print 'limit:', i, t
            limit = i
            break
    return limit
    
def cyto_export(tested, fname, negmult=100):
    ut.write_tab_file([(t[0],t[1],myml.rescale(t[2],negmult)) for t in
        tested], fname)
