from __future__ import division
import numpy as np
import utils as ut
import pairdict as pd
import corum as co

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

def jaccard(a, b):
    """
    Jaccard overlap of sets a,b.
    """
    return len(set.intersection(a,b))/len(set.union(a,b))

def ppv(a,b):
    """ 
    I had called this asym_overlap.  Brohee 2006 calls it pos predictive value.
    """
    return len(set.intersection(a,b))/len(a)

def simpson(a,b):
    return len(set.intersection(a,b))/min([len(a),len(b)])

def ppv_adj(a,b):
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


        
