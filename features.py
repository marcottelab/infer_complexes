from __future__ import division
import multiprocessing
import numpy as np
import operator
import os
import random
import sys
import pairdict as pd
import orth
import utils as ut

NCORES = multiprocessing.cpu_count()
NET_SPS = 'HS CE DM SC'.split()
FIVE_SPS_PLUS_NET = "CE DM HS MM SP".split()

def boot_arr(arr):
    return arr[ut.sample_wr(range(len(arr)), len(arr))]

def match_cols(arr, colstring):
    """
    colstring like 'Hs Ce HS ...' or 'Hs Ce HS-GN ...'
    """
    return [f for f in arr.dtype.names if (f[:2] in colstring.split()
                                           or f in colstring.split())]

def unmatch_cols(arr, colstring):
    """
    colstring like 'Hs Ce HS ...' or 'Hs Ce HS-GN ...'
    """
    return [f for f in arr.dtype.names[3:] if (f[:2] not in colstring.split()
                                           and f not in colstring.split())]

def keep_rows(arr, rows):
    return arr_copy(arr[rows])

def keep_cols(arr, colnames):
    return arr_copy(arr[['id1','id2','hit'] + colnames])

def balance_train(arrfeats):
    pos_inds, neg_inds = [[i for i,r in enumerate(arrfeats) if r[2]==pn] 
            for pn in 1,0]
    use_negs = min(len(pos_inds), len(neg_inds))
    return arrfeats[pos_inds + random.sample(neg_inds, use_negs)]

def regex_cols(arr, pattern):
    return arr[['id1','id2','hit'] + ut.regex_filter(arr.dtype.names, pattern)]

def merge_features(arr, pattern, func, do_remove):
    feats = ut.regex_filter(arr.dtype.names, pattern)
    print "features matched:", feats
    if len(feats) > 1:
        merged = [func(tuple(i)) for i in arr[feats]]
        name = pattern + '_' + func.__name__
        removefeats = set(feats) if do_remove else set([])
        keepdtypes = [dt for dt in arr.dtype.descr if not dt[0] in removefeats]
        newdtype = np.dtype(keepdtypes + [(name, arr[feats[0]].dtype)])
        newarr = np.empty(arr.shape, dtype=newdtype)
        for field in [d[0] for d in keepdtypes]:
            newarr[field] = arr[field]
        newarr[name] = merged
        return newarr
    else:
        print "Not enough features to merge."
        return arr

def retype_arr(arr, oldtype='f2', newtype='f4'):
    newdtype = [(dt[0],dt[1].replace(oldtype,newtype)) 
            for dt in arr.dtype.descr]
    return arr_copy(arr, newdtype=newdtype)

def arr_copy(arr, newdtype=None):
    dtype = newdtype if newdtype else arr.dtype
    newarr = np.empty(arr.shape, dtype=dtype)
    for field in newarr.dtype.names:
        newarr[field] = arr[field]
    return newarr
    
def arr_kfold(arr, kfold, k):
    total = arr.shape[0]
    start,stop = int(k*total/kfold), int((k+1)*total/kfold)
    print 'k: %s   start: %s   stop: %s' % (k, start, stop)
    arrtrain = keep_rows(arr, range(0,start)+range(stop,arr.shape[0]))
    arrtest = keep_rows(arr, range(start,stop))
    return arrtrain,arrtest

def merge_by_species(arr, matches, func, remove=False):
    """
    matches: like [apex] or [wcc, apex, ...]
    Makes patterns with match for each species, like 'Hs.*apex
    """
    assert not isinstance(matches,str), "matches is list, not string"
    def merge_recurse(arr, patterns, func):
        if patterns:
            newarr = merge_features(arr, patterns[0], func, remove)
            return merge_recurse(newarr, patterns[1:], func)
        else:
            return arr
    # Won't match a merged feature since that will have a * in it.
    patterns = [sp+'\w*'+match for match in matches for sp in
            ut.config()['elut_species'].split('_')]
    return merge_recurse(arr, patterns, func)

def species_cols(arr, set_species, count_ext=False):
    features = arr.dtype.names[3:]
    if set_species:
        cols = [f for f in features if f[:2] in set_species or (count_ext==True
            and f[:2]=='ex' and f[4:6] in set_species) ]
    else:
        cols = list(arr.dtype.names[3:])
    return cols

def filter_require_sp(arr, set_species, cutoff=0.5, count_ext=False):
    """
    Set_species: if None, just requires any column in the array to pass the
    cutoff.
    """
    cols = species_cols(arr, set_species, count_ext=count_ext)
    sp_max = [max(r) for r in arr[cols]]
    passing_inds = [i for i,m in enumerate(sp_max) if m > cutoff]
    return arr[passing_inds]

def filter_multi_orths(arr_in, basesp, cutoff):
    """
    For every interaction without base basespecies evidence, remove other basespecies
    evidence for that interaction when the base basespecies side of the orthogroup
    is greater than 1.
    """
    print "Filtering: require rows w/o %s > %s to have single orths" % (basesp,
            cutoff)
    arr = ut.arr_copy(arr_in)
    basesp_cols = [n for n in arr.dtype.names[3:] if n[:2]==basesp]
    assert len(basesp_cols)>0, 'No base species data.'
    maxes = arr_collist_maxes(arr, [basesp_cols])
    othersps = species_list(arr.dtype.names[3:])
    othersps.remove(basesp)
    spcols = [(sp, [n for n in arr.dtype.names[3:] if n[:2]==sp])
            for sp in othersps]
    ogs_all = orth.all_ogroup_sizes(basesp, othersps)
    cleared = 0
    for i in range(len(arr)):
        if maxes[i] < cutoff:
            row = arr[i]
            id1,id2 = row['id1'],row['id2']
            for sp, cols in spcols:
                ogsize_sp = ogs_all[sp]
                if (id1 in ogsize_sp and ogsize_sp[id1]>1) or (id2 in ogsize_sp
                        and ogsize_sp[id2]>1):
                    for col in cols: arr[i][col] = 0
                    cleared += 1
    print "%s species-sections of rows cleared" % cleared
    return arr

def maxes_fracs(arr):
    sps = species_list(arr.dtype.names[3:])
    fracs = set(['_'.join(n.split('_')[:-1]) for n in arr.dtype.names 
        if n[:2] in sps])
    fraclists = [[n for n in arr.dtype.names 
        if '_'.join(n.split('_')[:-1]) == frac] for frac in fracs]
    return arr_collist_maxes(arr, fraclists)

def arr_collist_maxes(arr, colsets):
    """
    Colsets is a list of lists of columns.  For each list of columns within
    colsets, we return a column in the new array that's the max of those.
    """
    maxes = np.zeros((len(arr), len(colsets)))
    for col in range(maxes.shape[1]):
        # Otherwise throws cryptic error if there's only one column
        assert min([len(cols)>1 for cols in colsets]) == True, "Single columns don't work."
        arrcols = arr[colsets[col]]
        maxfunc = max if len(colsets[col])>1 else lambda x: x
        maxes[:,col] = [maxfunc(arrcols[i]) for i in range(maxes.shape[0])]
    return maxes

def passing_scores(arr, thresh, counts=range(1,11), hits=None):
    hits = arr['hit'] if hits is None else hits
    npassing= [len([m for m in arr[i] if m>thresh]) for i in range(len(arr))]
    num_neg_pos = [len([i for i in hits if i==pn]) for pn in 0,1]
    proportions_neg_pos = [[len([i for i in range(len(arr)) if
        npassing[i]>=count and hits[i]==pn])/num_neg_pos[pn] for pn
        in 0,1] for count in counts]
    for i,(neg,pos) in zip(counts, proportions_neg_pos):
        print "%s: %0.2f p; %0.2f n; ratio %0.2f" % (i,pos,neg,pos/neg)
    return proportions_neg_pos

def max_fracs_passing(arr, thresh):
    maxes = maxes_fracs(arr)
    npassing= [len([m for m in maxes[i] if m>thresh]) 
            for i in range(len(maxes))]
    return npassing

def add_passing_features(arr, threshes):
    feat_names = ['passing%s' % t for t in threshes]
    newarr = ut.arr_add_feats(arr, feat_names)
    for t,name in zip(threshes, feat_names):
        newarr[name] = max_fracs_passing(arr, t)
    return newarr

def species_list(features, count_ext=False):
    sps = set([n[:2] for n in features 
                if n[:2] not in set(NET_SPS) and n[:2] != 'ex'])
    if count_ext and len(sps)==0:
        sps = set(FIVE_SPS_PLUS_NET)
    return sps

def passing_species(arr, cutoff=0.5, count_ext=False):
    """
    Intended for a single row of the array: a length=1 ndarray.
    """
    features = arr.dtype.names[3:]
    sps = species_list(features)
    spcols = [[f for f in features 
                if f[:2]==s or (count_ext==True and f[:2]=='ex' and f[4:6]==s)]
                for s in sps]
    maxes = arr_collist_maxes(arr, spcols)
    passing = [s for s,m in zip(sps, maxes[0]) if m > cutoff]
    return passing

def passing_fractionations_single(arr, cutoff=0.5):
    """
    Intended for a single row of the array: a length=1 ndarray. 
    """
    def frac_from_feat(f):
        f = f.replace('pq_euc','pq-euc')
        return '_'.join(f.split('_')[:-1]) 
    features = arr.dtype.names[3:]
    sps = species_list(features)
    fractionations = set([frac_from_feat(f) for f in features if f[:2] in sps])
    fcols = [[f for f in features if frac_from_feat(f)==frac] 
            for frac in fractionations]
    maxes = arr_collist_maxes(arr, fcols)
    #passing = [s for s,m in zip(sps, maxes[0]) if m > cutoff]
    return len([m for m in maxes[0] if m > cutoff])

def passing_fractionations(arr, cutoff=0.5):
    return [passing_fractionations_single(arr[i:i+1], cutoff=cutoff) 
            for i in range(len(arr))]

def passing_species_separate(arr, cutoff=0.5):
    """
    Intended for a single row of the array: a length=1 ndarray.
    """
    features = arr.dtype.names[3:]
    sps = species_list(features)
    spcols = [[f for f in features if f[:2]==s] for s in sps]
    spcols_ext = [[f for f in features if f[:2]=='ex' and f[4:6]==s]
                for s in sps]
    spcols_net = [[f for f in features if f[:2]==s.upper()] for s in sps]
    maxes = maxes_singles(arr, spcols)
    maxes_ext = maxes_singles(arr, spcols_ext)
    maxes_net = maxes_singles(arr, spcols_net)
    passing = ['frac' if m > cutoff else 'ext' if mext > cutoff else 'net' if
            mnet > cutoff else 'none'
        for s,m,mext,mnet in zip(sps, maxes, maxes_ext, maxes_net)]
    return passing

def arr_add_spsummary(arrfeats, cutoff):
    sps = species_list(arrfeats.dtype.names[3:])
    names = [s + '_evidence' for s in sps]
    newarr = ut.arr_add_feats(arrfeats, names, newtype='|S20')
    for i in range(len(newarr)):
        splist = passing_species_separate(arrfeats[i:i+1], cutoff)
        ut.arr_assign_row_values(newarr,i,range(len(arrfeats.dtype.names),
            len(newarr.dtype.names)), splist)
    return newarr

def arrfeats_add_ppis(arrfeats, ppis, name='ppi_score'):
    #assert len(arrfeats) == len(ppis), "Mismatched lengths"
    pdppis = pd.PairDict(ppis)
    newarr = ut.arr_add_feats(arrfeats, [name], newtype='<f2')
    notfound = 0
    for row in newarr:
        ppi = pdppis.find((row[0],row[1]))
        if ppi:
            row[-1] = pdppis.d[ppi][0]
        else:
            notfound += 1
    print "Not found:", notfound
    return newarr

def arrfeats_set_gold(arrfeats, pdgold):
    arrfeats = ut.arr_copy(arrfeats)
    for row in arrfeats:
        if pdgold.contains((row[0],row[1])):
            row[2] = 1
        else:
            row[2] = 0
    return arrfeats

def maxes_singles(arr, colsets):
    return [0 if len(cols)==0 else max(arr[cols][0]) for cols in colsets]

def filter_nsp(arr, nsp=2, cutoff=0.5, count_ext=False, do_counts=True):
    """
    Filter to only leave interactions for which evidence exists in nsp species.
    scores: [arr_train, arr_test]
    In addition to returning the filtered array, returns a matching pairdict of
    species counts.
    """
    features = arr.dtype.names[3:]
    sps = species_list(features, count_ext=count_ext)
    print 'Filtering >=%s of these species >=%s:' % (nsp, cutoff), sps
    spcols = [[f for f in features if f[:2]==s] for s in sps]
    if count_ext:
        spcols = [[f for f in features 
                if f[:2].lower()==s.lower() or (count_ext==True and f[:2]=='ex'
                    and f[4:6].lower()==s.lower())]
                for s in sps]
        spcols = [cols for cols in spcols if len(cols)>1]
        print "actually using these cols", spcols
    print '' if count_ext else 'NOT', 'Counting external in filter passing.'
    maxes = arr_collist_maxes(arr, spcols)
    exceed_inds = [i for i in range(len(arr))
                   if len([1 for m in maxes[i] if m > cutoff]) >= nsp]
    arrfilt = arr[exceed_inds]
    if do_counts:
        orig_sp_counts = [len([1 for m in maxes[i] if m > cutoff]) for i in
                range(len(arr))]
        sp_counts_filt = [c for c in orig_sp_counts if c >=nsp]
        pd_spcounts = pd.PairDict([[arrfilt[i][0],arrfilt[i][1],
            sp_counts_filt[i]] 
            for i in range(len(arrfilt))])
    else:
        pd_spcounts = None
    return arrfilt, pd_spcounts

def filter_nsp_nocounts(arr, **kwargs):
    return filter_nsp(arr, do_counts=False, **kwargs)[0]

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 6:
        sys.exit("usage: python myml.py f_examples f_classifier \
               perslice islice path")
    exs = np.load(sys.argv[1])
    clf = ut.loadpy(sys.argv[2])
    perslice = int(sys.argv[3])
    i = int(sys.argv[4])
    path = sys.argv[5]
    exs_slice = exs[i*perslice:(i+1)*perslice]
    del exs
    classify_slice(clf, exs_slice, 100000, savef=(path+str(i)+'_'),
            maintain=False, startslice=0)
