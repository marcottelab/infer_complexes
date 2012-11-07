from __future__ import division
import numpy as np
import utils as ut
import pairdict as pd
import corum as co
import cv
import operator
from Struct import Struct
import itertools as it


cp_funcs = [
    # 0 count of partially recovered gold cxs, first measure of havig/hart
    ("bader", lambda gold,cxs: count_overlaps(gold, cxs, limit=0.25,
        func=bader_score)),
    # 1 clustering-wise sensitivity, Brohee 2006
    ("sensitivity", lambda gold,cxs: avg_best(cxs, gold,
        sensitivity)),
    # 2 clustering-wise ppv, Brohee 2006
    ('ppv', lambda gold,cxs: ppv(gold, cxs)),
    # 3 geometric accuracy, Brohee 2006, third havig/hart
    ('accuracy', lambda gold,cxs: geom_avg(avg_best(cxs, gold,
        sensitivity), ppv(gold, cxs))),
    # 4 geometric separation, Brohee 2006
    ('separation', lambda gold, cxs: geom_avg(separation(gold, cxs),
        separation(cxs, gold))),
    ('overlaps sens', lambda gold,cxs: (count_overlaps(gold, cxs, limit=0.5,
        func=sensitivity)/len(gold))),
    ('overlaps sens rev', lambda gold,cxs: (count_overlaps(cxs, gold, limit=0.5,
        func=sensitivity)/len(cxs))),
    ('corr', lambda gold,cxs: gold_corr(gold, cxs)),
    #lambda gold,cxs: avg_best(cxs, gold, sensitivity_adj)
    # 8 Nepusz 2012, second and optimized metric for havig/hart
    ('mmr',lambda gold,cxs: max_match_ratio(gold, cxs)),
    ]

def stats(gold,cxs_list, conv_to_sets=True, func_inds=None):
    if conv_to_sets:
        gold = [set(c) for c in gold]
        cxs_list = [[set(c) for c in cxs] for cxs in cxs_list]
    use_funcs = np.array(cp_funcs)[func_inds] if func_inds else cp_funcs
    names = ut.i0(use_funcs)
    print names
    arr = np.zeros(len(cxs_list),dtype=','.join(['f8']*len(names)))
    arr.dtype.names = names
    for i, cxs in enumerate(cxs_list):
        print 'complex %s of %s' %(i, len(cxs_list))
        arr[i] = np.array([f(gold, cxs) for f in ut.i1(use_funcs)])
    return arr

def result_stats(sp, splits, clusts, nsp, func_inds=[2,8], split_inds=[2]):
    clstats = stats(result_gold(splits, sp, split_inds=split_inds,
        consv_sp=('Dm' if nsp==2 else '')), [c[1][0] for c in clusts],
        func_inds=func_inds)
    return Struct(clusts=clusts, stats=clstats)


def select_best(clstruct, scorenames, rfunc=operator.add, use_norm=True):
    clusts, stats = clstruct.clusts, clstruct.stats
    if use_norm: stats = norm_columns(stats)
    inds = np.argsort(reduce(rfunc, [stats[n] for n in scorenames]))[::-1]
    for i in inds[:10]: 
        print i, ["%0.4f " % s for s in clstruct.stats[i]], len(clusts[i][1][0]), len(clusts[i][1][1]), clusts[i][0]
    return clusts[inds[0]][1][0], clusts[inds[0]][1][1], inds[0]

def result_gold(splits, species, split_inds, make_unmerged=False,
        consv_sp='Dm'):
    if make_unmerged: 
        print "Converting to unmerged using conserved:", (consv_sp if consv_sp
                else "None")
        ppi_corum,_,_ = ppi.load_training_complexes(species,consv_sp)
        splits = unmerged_splits_from_merged_splits(ut.i1(ppi_corum), splits)
    gold = ut.i1(reduce(operator.add, [splits[i] for i in split_inds]))
    return gold

def prstats(gold, pred_ints, ntest_pos, prec=0.5):
    """
    Pred_ints: list of pred_result.ppis, or clustering filtered predicted ppis.
    """
    gold_pd = pd.PairDict([(p0,p1,1) for p0,p1 in
        co.pairs_from_complexes(dict([(i,g) for i,g in enumerate(gold)]))])
    tested_ppis = [[(id1,id2,score,1 if gold_pd.contains((id1,id2)) else
        0) for id1,id2,score,_ in ppis] for ppis in pred_ints]
    stats = [cv.preccheck(c,pchecks=[prec],total_trues=ntest_pos)[0]
            for c in tested_ppis]
    return stats

def unmerged_splits_from_merged_splits(unmerged, merged_splits):
    """
    Unmerged: list of sets.
    Merged_splits: A list of sets for each split (often 3 splits).
    """
    usplits = [[] for ms in merged_splits]
    for u in unmerged:
        which_split = np.argmax([best_match(u, ms, sensitivity) 
            for ms in merged_splits])
        usplits[which_split].append(u)
    return usplits

def norm_columns(arr):
    newarr = ut.arr_copy(arr)
    for n in newarr.dtype.names:
        newarr[n] = newarr[n]/np.max(newarr[n])
    return newarr

def max_match_ratio(gold, cxs):
    """
    Nepusz 2012, deciding score used for c1 param choice in havig/hart.
    """
    arr = arr_intersects(gold, cxs, func=bader_score)
    inds = zip(*np.where(arr>0)) # pairs of complexes that are matching at all.
    ind_scores = [(ind,arr[ind[0],ind[1]]) for ind in inds]
    ind_scores.sort(key=lambda x:x[1], reverse=True)
    x_used, y_used = set([]),set([])
    scores = []
    for (x,y),score in ind_scores:
        if x not in x_used and y not in y_used:
            x_used.add(x); y_used.add(y)
            scores.append(score)
    return sum(scores)/len(gold)


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


def arr_intersects(setsa, setsb, func=lambda a,b:len(set.intersection(a,b))):
    arr = np.zeros((len(setsa),len(setsb)))
    for i,a in enumerate(setsa):
        for j,b in enumerate(setsb):
            arr[i,j] = func(a,b)
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

def best_match_index(a, setsb, func):
    return np.argmax([func(a, complexb) for complexb in setsb])

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
    return pds_overlap([pd.PairDict(ppis) for ppis in pairs_iterables])

def pds_overlap(pds):
    """
    pds: list of pairdicts
    """
    return len(reduce(pd.pd_intersect_avals, pds).d)

def pds_alloverlaps(named_pds):
    """
    input: [(name, pairdict), ...]
    """
    for num in range(2,len(named_pds)+1):
        for n_pd in it.combinations(named_pds, num):
            print ut.i0(n_pd), pds_overlap(ut.i1(n_pd))

def triple_venn(three_ppis, names=['a','b','c']):
    ppis_names = zip(three_ppis, names)
    print zip(names, [len(p) for p in three_ppis])
    trip = ints_overlap(three_ppis)
    print names, trip
    intersects = []
    for (a,namea),(b,nameb) in it.combinations(ppis_names,2):
        intersect = ints_overlap([a,b])
        print namea, nameb, "--", intersect
        print namea, nameb, "-only-", intersect - trip







def cxs_self_match_frac(cxs, limit=.6, func=bader_score):
    return len(overlaps(cxs,cxs,limit,func=func,skip_self=True))/len(cxs)

