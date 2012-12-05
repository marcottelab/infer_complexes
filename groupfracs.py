from __future__ import division
from collections import defaultdict
import numpy as np
import scipy
import itertools as it
from Struct import Struct
import utils as ut
import elution as el
import score as sc

def load_normeluts(elutfs):
    return [el.NormElut(f) for f in elutfs]

class SplitScores(Struct):

    def __init__(self, split, list_subsplit_scores):
        self.split = split
        self.subscores = list_subsplit_scores

class SubsplitScore(Struct):

    def __init__(self, subsplit, cx_dist, excluded_dist, frac_name):
        self.subsplit = subsplit # split members actually observed
        self.cx_dist = cx_dist
        self.excluded_dist = excluded_dist
        self.frac = frac_name

def summarize(gids, normed_eluts, id2name=None, cx_max=.8, apart_min=.9,
        sort_sp=None):
    whole = score_whole(gids, normed_eluts)
    whole_spcounts = defaultdict(int)
    for score in whole.subscores:
        sp = score.frac[:2]
        if score.cx_dist < cx_max:
            whole_spcounts[sp] += 1
    print 'whole: %s total fractionations; passing:' %len(whole.split), \
            whole_spcounts.items()
    splits_scores = score_splits(gids, normed_eluts, cx_max=cx_max,
            apart_min=apart_min)
    splits_counts = []
    for split_scores in splits_scores:
        split_spcounts = defaultdict(int)
        for score in split_scores.subscores:
            sp = score.frac[:2]
            split_spcounts[sp] += 1
        splits_counts.append((split_scores.split,split_spcounts))
    if sort_sp is not None:
        sortfunc = lambda x: x[1][sort_sp] 
    else:
        sortfunc = lambda x: sum([v for k,v in x[1].items()])
    splits_counts.sort(key = sortfunc, reverse=True)
    splits_counts.sort(key = lambda x: len(x[0][0]), reverse=True)
    for x in splits_counts[:10]:
        if id2name:
            print ut.lol_trans(x[0],id2name), x[1].items()
        else:
            print x[0], x[1].items()
    return whole_spcounts, splits_counts

def score_whole(gids, norm_eluts):
    scores = []
    for e in norm_eluts:        
        if min([i in e.baseid2inds for i in gids])==True:
            rows = [ind for gid in gids if gid in e.baseid2inds 
                for ind in e.baseid2inds[gid]] 
            score = SubsplitScore(gids, distance_set(e.normarr, rows), None,
                    ut.shortname(e.filename)) 
            scores.append(score)
    return SplitScores(gids, scores)

def score_splits(gids, normed_eluts, cx_max=.8, apart_min=.9):
    splits = make_splits(gids)
    scored = []
    for split in splits:
        subscores = score_eluts_in_out(normed_eluts, split)
        subscores = [sc for sc in subscores 
                if sc.cx_dist < cx_max and sc.excluded_dist > apart_min]
        if len(subscores)>0:
            scored.append(SplitScores(split, subscores))
    return scored

def score_eluts_in_out(normed_eluts, gid_split, sp='Hs', **kwargs):
    scores = []
    for e in normed_eluts:
        rowid_split = [[ind for gid in gids if gid in e.baseid2inds 
            for ind in e.baseid2inds[gid]] 
            for gids in gid_split]
        if len(rowid_split[0])>=2 and len(rowid_split[1])>=1:#skip if any empties
            cx_dist = distance_set(e.normarr, rowid_split[0])
            exclude_dist = distance_exclude(e.normarr, rowid_split[0],
                    rowid_split[1])
            used_ids = [[e.prots[r] for r in rows] for rows in rowid_split]
            subscore = SubsplitScore(used_ids, cx_dist, exclude_dist,
                    ut.shortname(e.filename))
            scores.append(subscore)
    return scores

def distance_exclude(arr, ingroup, outgroup, combfunc=np.mean):
    """
    Ideally every member of outgroup is distant from every member of ingroup.
    """
    return combfunc([distance_pair(arr, i, j) for i in ingroup for j in outgroup])

def max_row(arr):
    maxind = np.argmax([np.sum(arr[i]) for i in range(arr.shape[0])])
    return arr[maxind]

def distance_set(arr, inds, combfunc = max_row):
    dists = np.zeros((scipy.misc.comb(len(inds),2, exact=1), arr.shape[1]))
    for pairind, (i,j) in enumerate(it.combinations(inds,2)):
        dists[pairind] = np.abs(arr[i]-arr[j])
    #print np.round(dists,2)
    comb_dists = combfunc(dists)
    #print np.round(comb_dists,2)
    dist = np.sum(comb_dists) / 2 # 2: Assumes rows are all normalized to 1.
    #print np.round(dist,2)
    return dist

def distance_pair(arr, ind1, ind2, **kwargs):
    x,y = arr[ind1], arr[ind2]
    dist = np.sum(np.abs(x-y))/2 # 2: Assumes rows are normalized to 1
    #print dist, np.sum(x), np.sum(y)
    return dist

def distance_group(arr, inds):
    return np.sum([distance_pair(arr, i, j) for i,j in it.combinations(inds,2)])

def max_cols(arr):
    return np.max(arr, axis=0)


def make_splits(items):
    set_all = set(items)
    # size 1 complexes are meaningless--start from 2
    sets_in = [set(x) for i in range(2,len(items)) 
            for x in it.combinations(items,i)]
    splits = [(s, set_all-s) for s in sets_in]
    return splits

def cos_score(arr, inds, **kwargs):
    return sum(reduce(np.multiply, [arr[i] for i in inds]))

def min_mult_score(arr, inds, power=2, **kwargs):
    return sum(np.min(arr[inds],axis=0)) * len(inds)**power

def score_eluts_dep(normed_eluts, gids, score_func=min_mult_score, **kwargs):
    score = 0
    for e in normed_eluts:
        if min([gid in e.pinv for gid in gids]) == True:
            score += score_func(e.normarr, [e.pinv[gid] for gid in gids],
                    **kwargs)
    return score

def r2(num):
    return round(num,2)
