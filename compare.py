from __future__ import division
import utils as ut
import pairdict as pd

def jaccard(a, b):
    """
    Jaccard overlap of sets a,b.
    """
    return len(set.intersection(a,b))/len(set.union(a,b))

def asym_overlap(a,b):
    """ 
    Related but not quite Simpson's Coefficient--would use min(len(a),len(b))
    as the denominator.
    """
    return len(set.intersection(a,b))/len(a)

def first_match(a, setsb, limit, skip_self):
    if skip_self:
        newsetsb = list(setsb)
        newsetsb.remove(a)
        setsb = newsetsb
    for b in setsb:
        if asym_overlap(a,b) > limit:
            return b

def matches(a, setsb, limit=0.5):
    for b in setsb:
        if jaccard(a,b) > limit:
            yield b

def overlaps(setsa, setsb, limit, skip_self):
    #return [a for a in setsa if matches(a,setsb).next()]
    return [a for a in setsa if first_match(a,setsb, limit, skip_self)]

def count_overlaps(setsa, setsb, limit=0.5, skip_self=False):
    return len(list(overlaps(setsa,setsb, limit,skip_self)))

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

