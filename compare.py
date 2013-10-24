from __future__ import division
import numpy as np
import scipy
import utils as ut
import pairdict as pd
import ppi
import corum as co
import cv
import operator
from Struct import Struct
import itertools as it
import sys
from pandas import DataFrame


cp_funcs = [
    # ( function name, function, whether to use random mean/variance )
    # 0 count of partially recovered gold cxs, first measure of havig/hart
    ("bader", lambda gold,cxs,cxppis: count_overlaps(gold, cxs, limit=0.25,
        func=bader_score), False),
    # 1 clustering-wise sensitivity, Brohee 2006
    ("sensitivity", lambda gold,cxs,cxppis: avg_best(cxs, gold,
        sensitivity), True),
    # 2 clustering-wise ppv, Brohee 2006
    ('ppv', lambda gold,cxs,cxppis: ppv(gold, cxs), True),
    # 3 geometric accuracy, Brohee 2006, third havig/hart
    ('accuracy', lambda gold,cxs,cxppis: geom_avg(avg_best(cxs, gold,
        sensitivity), ppv(gold, cxs)), True),
    # 4 geometric separation, Brohee 2006
    ('separation', lambda gold,cxs,cxppis: geom_avg(separation(gold, cxs),
        separation(cxs, gold)), True),
    ('overlaps sens', lambda gold,cxs,cxppis: (count_overlaps(gold, cxs, limit=0.5,
        func=sensitivity)/len(gold)), True),
    ('overlaps sens rev', lambda gold,cxs,cxppis: (count_overlaps(cxs, gold, limit=0.5,
        func=sensitivity)/len(cxs)) if len(cxs)!=0 else 0, True),
    ('corr', lambda gold,cxs,cxppis: gold_corr(gold, cxs), True),
    #lambda gold,cxs,cxppis: avg_best(cxs, gold, sensitivity_adj)
    # 8 Nepusz 2012, second and optimized metric for havig/hart
    ('mmr',lambda gold,cxs,cxppis: max_match_ratio(gold, cxs), True),
    #('non_overlap', lambda gold,cxs,cxppis: 1-cxs_self_match_frac(cxs), False), False),
    #('auroc', lambda gold,cxs,cxppis: auroc(gold,cxppis), False), 
    ('aupr', lambda gold,cxs,cxppis: aupr(gold, cxppis), False), 
    ('nonov_iter', lambda gold,cxs,cxppis: 1-cxs_self_match_frac_iter(cxs), False), 
    ('cliqueness_3_20', lambda gold,cxs,cxppis: cliqueness(filter_len(cxs, 3, 20), cxppis), False),
    #('ppi_pr', lambda gold,cxs,cxppis: prstats(gold, cxppis), False), 
    ('sm_cxs_3', lambda gold,cxs,cxppis: 1/cxs_avg_size(filter_len(cxs, 3, None)), False),
    ('n_proteins', lambda gold,cxs,cxppis: count_uniques(cxs), False),
    ('n_complexes_3_20', lambda gold,cxs,cxppis: len(filter_len(cxs, 3, 20)), False)
]


def stats(gold,cxs_list, cxppis_list=None, conv_to_sets=True, funcs=None):
    if conv_to_sets:
        gold = [set(c) for c in gold]
        cxs_list = [[set(c) for c in cxs] for cxs in cxs_list]
    funcs = funcs or ut.i0(cp_funcs)
    use_funcs = [f for f in cp_funcs if f[0] in funcs]
    print funcs
    arr = np.zeros(len(cxs_list),dtype=','.join(['f8']*len(funcs)))
    arr.dtype.names = funcs
    print '%s total maps.' % len(cxs_list)
    for i, (cxs, cxppis) in enumerate(zip(cxs_list, cxppis_list)):
        #sys.stdout.write(str(i))
        print i
        arr[i] = np.array([f(gold, cxs, cxppis) for f in ut.i1(use_funcs)])
    return arr

def result_stats(sp, splits, cxstructs, nsp, funcs=None,
        split_inds=[0], cxs_sets=None, min_gold_size=3):
    """
    split_inds: specifies which of the cxs_splits to use as gold--merges these
    indices. mar 2013 there are just two splits now, the train/cv set and the
    holdout set, so use index 0 for selection and 1 for testing. Note tree
    overfits and performs especially well on the training set.
    """
    cxs_gold = result_gold(splits, sp, split_inds=split_inds, 
            consv_sp=('Dm' if nsp==2 else ''))
    if min_gold_size>1:
        cxs_gold = [c for c in cxs_gold if len(c) >= min_gold_size]
    clstats = stats(cxs_gold, cxs_sets or [cstr.cxs for cstr in
            cxstructs], [cstr.cxppis for cstr in cxstructs], funcs=funcs)
    return Struct(cxstructs=cxstructs, stats=clstats)


def select_best(clstruct,
scorenames=['sensitivity','mmr','aupr','cliqueness_3_20','nonov_iter','n_proteins','n_complexes_3_20'],
        rfunc=operator.add, use_norm=False, dispn=15, score_factors=None,
        use_ranks=True, output_ranks=False, print_ranks=False,
        require_scores=None):
    cxstructs, stats = clstruct.cxstructs, clstruct.stats
    clusts = [cxstr.cxs for cxstr in cxstructs]
    scorenames = scorenames or list(stats.dtype.names)
    stats = stats[scorenames]
    ranks = rank_columns(stats)
    if use_ranks:
        stats = ranks
    else:
        if use_norm: stats = norm_columns(stats)
        if score_factors: stats = rescale_columns(stats, score_factors)
    inds = np.argsort(reduce(rfunc, [stats[n] for n in scorenames]))[::-1]
    if require_scores is not None:
        for req_name,thresh in require_scores:
            thresh = (np.median(clstruct.stats[req_name]) if thresh is None
                    else thresh)
            inds = [i for i in inds if clstruct.stats[req_name][i] > thresh]
    nstats = len(stats)
    def filt_params(s):
        return " ".join([p[:2]+p.split('=')[1] for p in s.split(',')])
    show_columns = (scorenames if require_scores is None else
            scorenames+ut.i0(require_scores))
    d = DataFrame(clstruct.stats[inds[:dispn]][show_columns],
            index=["#%s: %sc %si %s" %
                (i,len(clusts[i]),len(cxstructs[i].cxppis),
                    filt_params(cxstructs[i].params)) for i in inds[:dispn]])
    print d.head(dispn)
    for i in inds[:dispn]: 
        #print (i, ["%0.4f " % s for s in clstruct.stats[i]], len(clusts[i]), 
                #len(cxstructs[i].cxppis), cxstructs[i].params)
        if print_ranks:
            print i, [nstats-s for s in ranks[i]]
    if output_ranks:
        return inds
    else:
        return clusts[inds[0]], cxstructs[inds[0]].cxppis, inds[0]

def result_gold(splits, species, split_inds, make_unmerged=False,
        consv_sp='Dm'):
    if make_unmerged: 
        print "Converting to unmerged using conserved:", (consv_sp if consv_sp
                else "None")
        ppi_corum,_,_ = ppi.load_training_complexes(species,'',consv_sp)
        splits = unmerged_splits_from_merged_splits(ut.i1(ppi_corum),
                [[ut.i1(s) for s in split] for split in splits])
    gold = ut.i1(reduce(operator.add, [splits[i] for i in split_inds]))
    return gold

def prstats(gold, ppis, prec=0.2):
    """
    Pred_ints: list of pred_result.ppis, or clustering filtered predicted ppis.
    """
    tested,ntest_pos = tested_ppis(gold, ppis)
    recalled = cv.preccheck(tested,pchecks=[prec],total_trues=ntest_pos)[0]
    return recalled

def aupr(gold, ppis):
    """
    Area under precision-recall.
    """
    tested,ntest_pos = tested_ppis(gold, ppis)
    return cv.aupr(tested, ntest_pos)

def unmerged_splits_from_merged_splits(unmerged, merged_splits):
    """
    Unmerged: list of enumerated/named complex sets, usually from
    ppi.load_training_complexes.
    Merged_splits: A list of sets for each split (often 3 splits).
    """
    merged_splits = [ut.i1(ms) for ms in merged_splits]
    unmerged = ut.i1(unmerged)
    usplits = [[] for ms in merged_splits] + [[]] # extra for non matches
    for u in unmerged:
        matches = [best_match(u, ms, sensitivity) for ms in merged_splits]
        which_split = np.argmax(matches) if max(matches) > 0 else -1
        usplits[which_split].append(u)
    return usplits

def gold_label_ppis(ppis, merged_splits, sp, gold_nsp):
    gold_consv = 'Dm' if gold_nsp>1 else ''
    ppi_cxs,_,_ = ppi.load_training_complexes(sp, '', gold_consv)
    train_cxs = unmerged_splits_from_merged_splits(ppi_cxs, merged_splits)[0]
    ppis = cv.gold_label_ppis(ppis, co.pairs_from_complexes(train_cxs))
    return ppis

def rescale_columns(arr, scale_factors):
    newarr = ut.arr_copy(arr)
    for i,n in enumerate(newarr.dtype.names):
        #newarr[n] = np.nan_to_num(newarr[n]/np.max(np.nan_to_num(newarr[n])))
        newarr[n] = newarr[n] * scale_factors[i]
    return newarr

def norm_columns(arr):
    newarr = ut.arr_copy(arr)
    for n in newarr.dtype.names:
        newarr[n] = scipy.stats.zscore(np.nan_to_num(newarr[n]))
    return newarr

#def rescale_rand(arr, randarr):
    #newarr = ut.arr_copy(arr)
    #for n in newarr.dtype.names:
        #randmean = np.mean(randarr[n])
        #datamean = np.mean(arr[n])
        #newarr[n] = ( newarr[n] 

def rank_columns(arr):
    newarr = ut.arr_copy(arr)
    for n in newarr.dtype.names:
        newarr[n] = scipy.stats.rankdata(np.nan_to_num(newarr[n]))
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
    ints_dict = co.corum_ints_duped([(i,set(ps)) 
        for i,ps in enumerate(cxs)])
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
    clust_ppvs = np.nan_to_num(np.max(arr, axis=0) / np.sum(arr, axis=0))
    clust_sizes = [len(c) for c in setsb]
    overall_ppv = weighted_avg(clust_ppvs, clust_sizes)
    #return np.sum(np.max(arr, axis=0))/np.sum(arr)
    return overall_ppv

def separation(setsa, setsb):
    # Row-wise separation. Brohee 2006.
    # Transpose arguments for col-wise separation.
    if len(setsa)==0 or len(setsb)==0: return 0
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

def subset_plus_one(a,b,minsize):
    if len(a)>=minsize and len(b)>=minsize:
        intersect = len(set.intersection(a,b))
        if intersect > 0:
            if min(len(a),len(b)) <= intersect + 1:
                return 1
    return 0

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

def avg_best(setsa, setsb, func, weighted=True):
    scores = [best_match(a, setsb, func) for a in setsa]
    if weighted:
        lengths = [len(c) for c in setsa]
        return weighted_avg(scores, lengths)
    else:
        return np.mean([scores])

def weighted_avg(scores, weights):
    return np.sum(np.array(scores)*np.array(weights)) / np.sum(weights)

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
    return len(pds_intersect_pd(pds).d)

def pds_intersect_pd(pds):
    """
    pds: list of pairdicts
    """
    return reduce(pd.pd_intersect_avals, pds)

def pds_alloverlaps(named_pds):
    """
    input: [(name, pairdict), ...]
    """
    for num in range(2,len(named_pds)+1):
        for n_pd in it.combinations(named_pds, num):
            print ut.i0(n_pd), pds_overlap(ut.i1(n_pd))

def triple_venn(three_ppis, names=['a','b','c']):
    # Can send output to venn3(sizes, names) for plotting
    ppis_names = zip(three_ppis, names)
    full_sizes = [len(p) for p in three_ppis]
    print zip(names, full_sizes)
    trip = ints_overlap(three_ppis)
    print names, trip
    intersects_2 = []
    for (a,namea),(b,nameb) in it.combinations(ppis_names,2):
        intersect = ints_overlap([a,b])
        print namea, nameb, "--", intersect
        intersect_2 = intersect - trip
        print namea, nameb, "-only-", intersect_2
        intersects_2.append((set([namea,nameb]), intersect_2))
    only_singles = []
    for i,(a,namea) in enumerate(ppis_names):
        #print namea, len(a), trip, intersects_2 #debug
        only_single = (len(a) - trip - sum([x[1] for x in intersects_2 if namea in x[0]]))
        print namea, "only:", only_single
        only_singles.append(only_single)
    # for output into matplotlib_venn format
    set_sizes = [0] * 7 
    set_sizes[6] = trip
    set_sizes[:2], set_sizes[3] = only_singles[:2], only_singles[2]
    set_sizes[2], set_sizes[4], set_sizes[5] = ut.i1(intersects_2)
    return set_sizes, names

def cxs_avg_size(cxs):
    return np.mean([len(c) for c in cxs])

def cxs_self_match_frac(cxs, limit=.6, func=bader_score):
    return (len(overlaps(cxs,cxs,limit,func=func,skip_self=True))/len(cxs) if
            len(cxs)!=0 else 0)

def cxs_self_match_frac_iter(cxs, ntests=10, func=bader_score):
    # should instead consider limit=1 rather than 0--1 is meaningful, 0 is not
    limits = np.arange(0,1,1/ntests)
    return np.mean([cxs_self_match_frac(cxs, limit=lim, func=func) for lim in
        limits])

def ints_overlap_orth(aints, bints, b2a):
    """
    Make a list of orthogroup pair sets, instead of just pairs of genes, then
    query those for a given pair of genes in species A.
    """
    def matchpairset(pair,losetpairs):
        for p in losetpairs:
            if ((pair[0] in p[0] and pair[1] in p[1]) or 
                    (pair[1] in p[0] and pair[0] in p[1])):
                return True
    bints_spa = [(b2a[g[0]],b2a[g[1]]) for g in bints 
            if g[0] in b2a and g[1] in b2a]
    return len([pair for pair in aints if (matchpairset(pair,bints_spa))])

def consv_pairs(ints, odict):
    return [p for p in ints if p[0] in odict and p[1] in odict]

def ints_overlap_consv(intlists, odict):
    """
    Get the overlap just of the interactions for which each member of the pair
    has an ortholog in the given odict.
    """
    return ints_overlap([consv_pairs(ints,odict) for ints in intlists])

def count_trues(ppis, gold, ncounts=10000):
    pdgold = pd.PairDict(gold)
    print "Total trues:", len(gold)
    count = 0
    for n,p in enumerate(ppis):
        if pdgold.contains((p[0],p[1])):
            count += 1
        if n % ncounts == 0:
            print "%s True of first %s" %(count, n)

def combine_clstructs(ca, cb):
    return Struct(cxstructs = ca.cxstructs + cb.cxstructs,
            stats=np.concatenate((ca.stats, cb.stats)))

def tested_ppis(gold_cxs, ppis):
    gold_ints = co.pairs_from_complexes(gold_cxs)
    ntest_pos = len(gold_ints)
    pdtrues = pd.PairDict(gold_ints)
    ppis = [(p[0],p[1],p[2],1 if pdtrues.contains(tuple(p[:2])) else 0) for p in
            ppis]
    return ppis, ntest_pos

def auroc(gold_cxs, cxppis):
    cxppis, ntest_pos = tested_ppis(gold_cxs, cxppis)
    xs,ys = cv.roc(cxppis) 
    return cv.auroc(xs,ys)

def cliqueness(cxs, cxppis):
    pdints = pd.PairDict(cxppis)
    return np.mean([clique_score(c,pdints) for c in cxs])

def clique_score(cx, pdints):
    cx_ints = co.pairs_from_complexes([cx])
    return len([1 for edge in cx_ints if pdints.contains(edge)])/len(cx_ints)

def filter_len(lol, min_len=0, max_len=None):
    return [items for items in lol if len(items) >= min_len and (max_len is None
        or len(items) <= max_len)]

def count_uniques(lol):
    return (len(reduce(set.union, [set(items) for items in lol])) if len(lol)>0
            else 0)
