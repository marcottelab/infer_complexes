from __future__ import division
import numpy as np
import utils as ut
import pairdict as pd
import corum as co


cp_funcs = [
    # count of partially recovered gold cxs, first measure of havig/hart
    lambda gold,cxs: (count_overlaps(gold, cxs, limit=0.25,
        func=bader_score)),
    # clustering-wise sensitivity, Brohee 2006
    lambda gold,cxs: avg_best(cxs, gold, sensitivity),
    # clustering-wise ppv, Brohee 2006
    lambda gold,cxs: ppv(gold, cxs),
    # geometric accuracy, Brohee 2006, second havig/hart
    lambda gold,cxs: geom_avg(avg_best(cxs, gold, sensitivity),
        ppv(gold, cxs)),
    # geometric separation, Brohee 2006
    lambda gold, cxs: geom_avg(separation(gold, cxs), separation(cxs,
        gold)),
    lambda gold,cxs: (count_overlaps(gold, cxs, limit=0.5,
        func=sensitivity)/len(gold)),
    lambda gold,cxs: (count_overlaps(cxs, gold, limit=0.5,
        func=sensitivity)/len(cxs)),
    lambda gold,cxs: gold_corr(gold, cxs),
    #lambda gold,cxs: avg_best(cxs, gold, sensitivity_adj)
    ]

def stats(gold,cxs_list, conv_to_sets=True, norm=True):
    if conv_to_sets:
        gold = [set(c) for c in gold]
        cxs_list = [[set(c) for c in cxs] for cxs in cxs_list]
    arr = np.zeros((len(cxs_list), len(cp_funcs)))
    for i, cxs in enumerate(cxs_list):
        arr[i,:] = [f(gold, cxs) for f in cp_funcs]
    if norm:
        arr = arr/np.max(arr,axis=0)
    return arr

def gold_corr(gold, cxs):
    """
    gold: List of complexes from the gold standard
    cxs: List of complexes to compare against gold standard
    """
    # Get list of all the proteins from the gold standard and make a nprots x
    # nprots array holding whether the proteins interact or not.  Fill this in
    # with the gold standard, and with the cxs to test.  Take the correlation
    # between these arrays.
    gold_prots = reduce(set.union, gold)
    gold_arr, cxs_arr = [complex_arr(c, gold_prots) for c in gold, cxs]
    return np.corrcoef(gold_arr.flat, cxs_arr.flat)[0,1]

def geom_avg(a,b):
    return np.math.sqrt(a*b)

def complex_arr(cxs, prots):
    arr = np.zeros((len(prots),len(prots)))
    ints_dict = co.corum_ints_duped(dict([(i,set(ps)) 
        for i,ps in enumerate(cxs)]))
    p_inds = ut.list_inv_to_dict(prots)
    for p,partners in ints_dict.items():
        if p in p_inds:
            for partner in partners:
                if partner in p_inds:
                    arr[p_inds[p], p_inds[partner]] = 1
    return arr

def bader_score(a,b):
    # Bader/Hogue 2003 from Havig/Hart 2012
    intersect = len(set.intersection(a,b))
    return intersect**2/(len(a)*len(b))

def jaccard(a, b):
    return len(set.intersection(a,b))/len(set.union(a,b))

def sensitivity(a,b):
    # Brohee 2006 uses this as the criteria for 'sensitivity'. I had called it
    # asym_overlap previously.
    return len(set.intersection(a,b))/len(a)

def ppv(setsa, setsb):
    # clustering-wise ppv from Brohee 2006
    arr = arr_intersects(setsa, setsb)
    # Redundant steps in how their algo is presented were elminated but kept
    # here for clarity
    #clust_ppvs = np.max(arr, axis=0) / np.sum(arr, axis=0)
    #overall_ppv = (np.sum(arr, axis=0)*clust_ppvs) / np.sum(arr)
    return np.sum(np.max(arr, axis=0))/np.sum(arr)

def separation(setsa, setsb):
    # Row-wise separation. Brohee 2006.
    # Transpose arguments for col-wise separation.
    arr = arr_intersects(setsa, setsb)
    col_normed = np.nan_to_num(arr/np.sum(arr, axis=0))
    row_normed = np.nan_to_num(np.array([arr[i,:]/np.sum(arr,axis=1)[i] 
        for i in range(arr.shape[0])]))
    comb = col_normed * row_normed
    return np.mean(np.sum(comb, axis=1))


def arr_intersects(setsa, setsb):
    arr = np.zeros((len(setsa),len(setsb)))
    for i,a in enumerate(setsa):
        for j,b in enumerate(setsb):
            arr[i,j] = len(set.intersection(a,b))
    return arr

def simpson(a,b):
    return len(set.intersection(a,b))/min([len(a),len(b)])

def sensitivity_adj(a,b):
    """ 
    Tries to remove the issue that smaller complexes get biased towards high
    scores by subtracting the base case where just a single member overlaps.
    """
    return (len(set.intersection(a,b))-1)/len(a)

def best_match(a, setsb, func):
    return max([func(a, complexb) for complexb in setsb])

def best_match_item(a, setsb, func):
    return max(setsb, key=lambda b: func(a,b))

def avg_best(setsa, setsb, func):
    return np.mean([best_match(a, setsb, func) for a in setsa])

def first_match(a, setsb, limit, skip_self=False, func=ppv):
    if skip_self:
        newsetsb = list(setsb)
        newsetsb.remove(a)
        setsb = newsetsb
    for b in setsb:
        if func(a,b) > limit:
            return b

def matches(a, setsb, limit=0.5):
    for b in setsb:
        if jaccard(a,b) > limit:
            yield b

def overlaps(setsa, setsb, limit, **kwargs):
    #return [a for a in setsa if matches(a,setsb).next()]
    return [a for a in setsa if first_match(a,setsb, limit, **kwargs)]

def count_overlaps(setsa, setsb, limit=0.5, **kwargs):
    return len(list(overlaps(setsa,setsb, limit,**kwargs)))

def matchpairset(pair,losetpairs):
    for p in losetpairs:
        if (pair[0] in p[0] and pair[1] in p[1]) or (pair[1] in p[0] and pair[0] in p[1]):
            return True

def ints_overlap(pairs_iterables):
    """ 
    Provide a list: [pairsa, pairsb(, pairsc)]
    """
    return pds_overlap(pairs_iterables)

def pds_overlap(pds):
    """
    pds: list of pairdicts
    """
    return len(reduce(pd.pd_intersect_avals, pds).d)


        
