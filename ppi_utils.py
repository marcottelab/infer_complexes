from __future__ import division
import itertools as it
import pairdict as pd
import utils as ut

def groups_to_pairs(lol_groups):
    raw_pairs = ut.flatten([[x for x in it.combinations(group,2)] 
        for group in lol_groups])
    deduped = dedupe(raw_pairs)
    return deduped

def load_ppis(fname):
    return [(p[0], p[1], float(p[2]), int(p[3])) 
            for p in ut.load_tab_file(fname)]

def dedupe(pairlist):
    return pd.PairDict(pairlist).d.keys()

def intersect(pairsa, pairsb):
    pda = pd.PairDict(pairsa)
    intersect = [p for p in pairsb if pda.contains(p[:2])]
    return dedupe(intersect)

def unique_items(pairs):
    return list(set(ut.i0(pairs) + ut.i1(pairs)))

