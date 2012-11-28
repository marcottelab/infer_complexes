from __future__ import division
from Struct import Struct
import numpy as np
import scipy
import itertools as it
import utils as ut
import elution as el


def prep_eluts(elutfs):
    eluts = (el.load_elution(f) for f in elutfs)
    return [Struct(pinv = ut.list_inv_to_dict(e.prots), 
        arr = normalize_fracs(e.mat)) for e in eluts]

def score_eluts(normed_eluts, gids, score_func=min_mult_score, **kwargs):
    score = 0
    for e in normed_eluts:
        if min([gid in e.pinv for gid in gids]) == True:
            score += score_func(e.arr, [e.pinv[gid] for gid in gids], **kwargs)
    return score

def cos_score(arr, inds, **kwargs):
    return sum(reduce(np.multiply, [arr[i] for i in inds]))

def min_mult_score(arr, inds, power=2, **kwargs):
    return sum(np.min(arr[inds],axis=0)) * len(inds)**power

def euc_score(arr, ind1, ind2, **kwargs):
    #  Overlap btw values a,b for a slice: max(ai,bi) - abs(ai-bi)
    # Relative fraction overlap: overlap*2 / (sum(a)+sum(b))
    # Only pairs so far
    x,y = arr[ind1], arr[ind2]
    dist = np.sum(np.abs(x-y))
    #print dist, np.sum(x), np.sum(y)
    return 1 - dist/(np.sum(x)+np.sum(y))

def euc_group(arr, inds):
    return np.sum([euc_score(arr, i, j) for i,j in it.combinations(inds,2)])

def max_cols(arr):
    return np.max(arr, axis=0)

def euc_set(arr, inds, combfunc = max_cols):
    dists = np.zeros((scipy.misc.comb(len(inds),2, exact=1), arr.shape[1]))
    for pairind, (i,j) in enumerate(it.combinations(inds,2)):
        dists[pairind] = np.abs(arr[i]-arr[j])
    #print dists
    comb_dists = combfunc(dists)
    #print comb_dists
    dist = np.sum(comb_dists)
    #print dist
    return 1 - dist/np.sum([np.sum(arr[i]) for i in inds])

def normalize_fracs(arr, norm_rows=True, norm_cols=True):
    if norm_cols:
        # Normalize columns first--seems correct for overall elution profile if
        # you're calculating correlation-type scores
        arr = ut.arr_norm(arr, 0)
    if norm_rows:
        arr = ut.arr_norm(arr, 1)
    return arr

