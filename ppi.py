import examples as ex
import corum as co
import seqs
import orth
import utils as ut
import elution as el
from Struct import Struct
import numpy as np
import fnet
import score
from numpy import zeros,ndarray
from pairdict import PairDict

def learning_examples(species, elut_fs, scores=['poisson','wcc','apex'],
                      extdata=['net_Hs19', 'ext_Dm_guru','ext_Hs_malo'], 
                      splits=[0,.33,.66,1], neg_ratios=[10,40],
                      ind_cycle=[0,-1], cutoff=0.25, base_tt=None,
                      pos_splits=None, test_negs=None):
    """
    Species: 'Hs', 'Ce', ...
    Provide base_tt as [arrtest,arrtrain] to use saved examples.
    Give test_negs as 'all' to use full set of proteins from fs elution data.
    """
    train,test = pd_from_tt(base_tt) if base_tt else \
                 base_splits(species, elut_fs, splits, neg_ratios,
                      ind_cycle, test_negs, pos_splits)
    ntest_pos = len([v for k,v in test.d.items() if v[0]==1])
    print 'total test positives:', ntest_pos
    atrain,atest = [new_score_array(pdict, scores, elut_fs, extdata)
                    for pdict in train,test]
    train_test = [score_and_filter(arr, scores, elut_fs, cutoff, species,
                      extdata) for arr in atrain,atest]
    print 'done.', stats(train_test)
    return train_test, ntest_pos

def pd_from_tt(train_test):
    return [PairDict([[p[0],p[1],p[2]] for p in t]) for t in train_test]

def base_splits(species, elut_fs, splits, neg_ratios, ind_cycle,
                test_negs, pos_splits=None):
    ppi_cxs,clean_cxs = load_training_complexes(species)
    elut_prots = el.all_prots(elut_fs, species) if test_negs=='all' else None
    train,test = ex.base_examples(ppi_cxs, clean_cxs, elut_prots, splits,
        nratio_train=neg_ratios[0], nratio_test=neg_ratios[1],
        pos_splits=pos_splits, ind_cycle=ind_cycle)
    return train,test

def predict_all(species, elut_fs, scores=['poisson','wcc','apex'],
              extdata=['net_Hs19', 'ext_Dm_guru','ext_Hs_malo'], cutoff=0.25):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    pairs = el.all_filtered_pairs(elut_fs, scores, cutoff, species)
    print len(pairs.d), 'total interactions passing cutoff'
    for k in pairs.d: pairs.d[k] = -1 #interaction marked as unknown
    arr = new_score_array(pairs, scores, elut_fs, extdata)
    del pairs #lots of memory
    return score_and_filter(arr, scores, elut_fs, cutoff, species,
                extdata, must_filter=False)

def score_and_filter(arr, scores, elut_fs, cutoff, species, 
                     extdata, must_filter=True):
    print 'Scoring.'
    score.score_array_multi(arr, species, elut_fs, scores, cutoff)
    if must_filter and cutoff != -1:
        print 'Filtering.'
        columns = range(3,len(arr[0]))
        # note that this double-filters often. only an efficiency issue.
        arr = filter_arr(arr, columns, cutoff)
    for d in extdata: fnet.score_arr(arr, species, d)
    return arr 

def new_score_array(pd, scores, fs, extdata):
    """
    Pre-allocate an array sized based on the elution files and extdata.
    arr[0] returns the first row of data.
    row[name_score(f,s)] returns the proper score
    row['id1','id2','hit'] index to the first three columns
    """
    data_names = [score.name_score(f,s) for f in fs for s in scores]
    for key in extdata:
        if key[:3] == 'net':
            data_names += fnet.fnet_names(ut.config()[key])
        else:
            data_names.append(key)
    arr = zeros(shape = len(pd.d),
                dtype = ','.join(['a20']*2 + ['i1'] + ['f2']*len(data_names)))
    names = ['id1','id2','hit'] + data_names
    arr.dtype.names = tuple(names)
    for i,((p1,p2),lhit) in enumerate(pd.d.items()):
        row = arr[i]
        row['id1'],row['id2'],row['hit'] = p1,p2,lhit
    return arr

def stats(train_test):
    nums =  [len([1 for row in arr if row['hit']==tf])
             for arr in train_test for tf in [1,0]]
    return 'train %sP/%sN; test %sP/%sN' % tuple(nums)

def load_training_complexes(species):
    ppi = co.load_complexes_singleline(ut.proj_path('ppi_cxs'))
    clean = co.load_complexes_multiline(ut.proj_path('clean_cxs'))
    ppi,clean = [convert_complexes(cxs, species) for cxs in [ppi,clean]]
    return ppi,clean

def convert_complexes(cxs, species):
    """
    Convert from corum uniprot to human default, then to other base species if
    specified.
    """
    dconv = ut.load_dict_sets(ut.proj_path('corum_convert'))
    ps = seqs.load_prots_from_fasta(ut.proj_path('fastadir', 'Hs.fasta'))
    cxs = co.convert_complexes(cxs, dconv, ps)
    if species!='Hs':
        dconv = orth.odict('Hs', species)
        ps = seqs.load_prots_from_fasta(ut.proj_path('fastadir', 
                                                    species+'.fasta'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    return cxs

def filter_arr(arr, columns, cutoff):
    """
    Uses max: return a new array of all rows from input arr with at least one
    value in columns exceeding cutoff.
    """
    newarr = arr[[row for row in range(len(arr)) if max([val for col,val in
                enumerate(arr[row]) if col in columns]) > cutoff]]
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

