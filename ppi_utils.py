from __future__ import division
import itertools as it
import random
import pairdict as pd
import utils as ut

def groups_to_pairs(lol_groups):
    raw_pairs = ut.flatten([[x for x in it.combinations(group,2)] 
        for group in lol_groups])
    deduped = dedupe(raw_pairs)
    return deduped

def load_ppis(fname):
    return ut.load_lol(fname, dtypes=(str,str,float,int))

def dedupe(pairlist):
    return pd.PairDict(pairlist).d.keys()

def intersect(pairsa, pairsb):
    pda = pd.PairDict(pairsa)
    intersect = [p for p in pairsb if pda.contains(p[:2])]
    return dedupe(intersect)

def unique_items(pairs):
    return set(ut.i0(pairs) + ut.i1(pairs))

def nonpairs_gen(pairs, n):
    items = list(set(ut.i0(pairs) + ut.i1(pairs)))
    exclude = set(pairs)
    pdexclude = pd.PairDict(exclude)
    count = 0
    while count < n:
        pair = (random.choice(items), random.choice(items))
        if not pdexclude.contains(pair):
            yield pair
            count += 1
