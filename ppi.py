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

def learning_examples(species, elut_fs, base_tt, nsp, scores=['poisson','wcc','apex'],
                      extdata=['net_Hs19', 'ext_Dm_guru','ext_Hs_malo'], 
                      splits=[0,.33,.66,1], neg_ratios=[2.5,230],
                      ind_cycle=[0,-1], cutoff=0.25, 
                      pos_splits=None, test_negs=None, gold_consv_sp='Dm', 
                      do_filter=True): 
    """
    - Species: 'Hs' or 'Ce'... The base of the predictions in terms of
      identifiers used and orthology pairings.
    - base_tt: as [arrtest,arrtrain] to use saved examples.
    - Test_negs: 'all' to use full set of proteins from fs elution data; None
      to generate test negatives like training negatives, from inter-complex
      interactions.
    - extdata: list of items that are either a string data_key matching to
      those specified in config.py, or a tuple ('name', data): made for _as_ext
      re-use of ml scores as a new feature, but could support any data set
      provided as a variable rather than file. The data should simply be a
      sequence of list/tuple of (id1, id2, score). other_evidence just appends
      to this.
    - do_filter: only set to false if you want a full unfiltered test/train arr
      """
    gold_consv_sp = gold_consv_sp if nsp>1 else '' #ignore for 1-sp case
    train,test = pd_from_tt(base_tt) if base_tt else \
                 base_splits(species, elut_fs, splits, neg_ratios,
                      ind_cycle, test_negs, pos_splits, gold_consv_sp)
    ntest_pos = len([v for k,v in test.d.items() if v[0]==1])
    # Note this is wrong if base_tt is supplied since that is already filtered.
    print 'total test positives:', ntest_pos
    atrain,atest = [new_score_array(pdict, scores, elut_fs, extdata) 
            for pdict in train,test]
    train_test = [score_and_filter(arr, scores, elut_fs, cutoff, species,
                      extdata, do_filter) for arr in atrain,atest]
    print 'done.', stats(train_test)
    return train_test, ntest_pos

def pd_from_tt(train_test):
    return [PairDict([[p[0],p[1],p[2]] for p in t]) for t in train_test]

def base_splits(species, elut_fs, splits, neg_ratios, ind_cycle,
                test_negs, pos_splits, consv_sp):
    ppi_cxs,clean_cxs = load_training_complexes(species, consv_sp)
    elut_prots = el.all_prots(elut_fs, species) if test_negs=='all' else None
    train,test = ex.base_examples(ppi_cxs, clean_cxs, elut_prots, splits,
        nratio_train=neg_ratios[0], nratio_test=neg_ratios[1],
        pos_splits=pos_splits, ind_cycle=ind_cycle)
    return train,test

def predict_all(species, elut_fs, scores=['poisson','wcc','apex'],
              extdata=['net_Hs19', 'ext_Dm_guru','ext_Hs_malo'], cutoff=0.25,
              base_arr=None, do_filter=False):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    if not base_arr:
        pd_all = el.all_filtered_pairs(elut_fs, scores, cutoff, species)
        print len(pd_all.d), 'total interactions passing cutoff'
        for k in pd_all.d: pd_all.d[k] = -1 #interaction marked as unknown
        arr = new_score_array(pd_all, scores, elut_fs, extdata)
        del pd_all #lots of memory
    else:
        arr = base_arr
    scored_arr = score_and_filter(arr, scores, elut_fs, cutoff, species,
                extdata, do_filter)
    return scored_arr

def score_and_filter(arr, scores, elut_fs, cutoff, species, extdata,
        do_filter):
    print '\nScoring %s elutions with %s base, scores: %s.' % (len(elut_fs),
            species, ','.join(scores))
    score.score_array_multi(arr, species, elut_fs, scores, cutoff)
    if cutoff != -1 and do_filter:
        print 'Filtering.'
        columns = range(3,len(arr[0]))
        # note that this double-filters often. only an efficiency issue.
        arr = filter_arr(arr, columns, cutoff)
    if extdata:
        print '\nScoring with external data:', ','.join(extdata)
        for d in extdata: fnet.score_arr_ext(arr, species, d)
    if cutoff == -1 and do_filter: # correct: filter this if we haven't before.
        arr = remove_empties(arr)
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
        ext_names = fnet.fnet_names(ut.config()[key])
        if ext_names:
            data_names += ext_names
        else:
            data_names.append(key)
    arr = base_array(((p1,p2,lhit) for (p1,p2),lhit in pd.d.items()),
            data_names, len(pd.d))
    return arr

def base_array(triple, data_names, data_len):
    """
    triple should be a sequence of (id1, id2, hit)
    """
    arr = zeros(shape = data_len, dtype = ','.join(['a20']*2 + ['i1'] +
        ['f2']*len(data_names)))
    names = ['id1','id2','hit'] + data_names
    arr.dtype.names = tuple(names)
    for i,(p1,p2,hit) in enumerate(triple):
        row = arr[i]
        row['id1'],row['id2'],row['hit'] = p1,p2,hit
    return arr

def stats(train_test):
    nums =  [len([1 for row in arr if row['hit']==tf])
             for arr in train_test for tf in [1,0]]
    return 'train %sP/%sN; test %sP/%sN' % tuple(nums)

def load_training_complexes(species, consv_sp):
    ppi = co.load_complexes_singleline(ut.proj_path('ppi_cxs'))
    clean = co.load_complexes_multiline(ut.proj_path('clean_cxs'))
    ppi,clean = [convert_complexes(cxs, species, consv_sp) for cxs in
            [ppi,clean]]
    return ppi,clean

def convert_complexes(cxs, species, consv_sp):
    """
    Convert from corum uniprot to human default, then to other base species if
    specified.
    If consv_sp (conserved species) is provided, restrict gold standard to only
    include interactions for proteins with orthologs in that species.
    """
    dconv = ut.load_dict_sets(ut.proj_path('corum_convert'))
    ps = seqs.load_prots_from_fasta(ut.proj_path('fastadir', 'Hs.fasta'))
    cxs = co.convert_complexes(cxs, dconv, ps)
    if species!='Hs':
        dconv = orth.odict('Hs', species)
        ps = seqs.load_prots_from_fasta(ut.proj_path('fastadir', 
                                                    species+'.fasta'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    if consv_sp:
        conserved = orth.odict(species, consv_sp)
        cxs = dict([(c,set([g for g in gs if g in conserved])) for c,gs in
            cxs.items()])
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

def remove_empties(arr):
    if len(arr)==2: # assume this is test-train
        return [remove_empties(a) for a in arr]
    newarr = arr[list(arr.dtype.names[3:])]
    return arr[np.where([max(newarr[i]) for i in range(len(newarr))])]
