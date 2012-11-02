from __future__ import division
import utils as ut
import pairdict as pd
import itertools as it

def poss_ints(nps):
    return (nps*(nps-1))/2

def cliqueness_mult(cxs,cxppis):
    pdppis = pd.PairDict(cxppis)
    return [cliqueness(c, pdppis) for c in cxs]

def cliqueness(cx,pdppis):
    return len([1 for p1,p2 in it.combinations(cx,2) if
        pdppis.contains((p1,p2))]) / poss_ints(len(cx)) 

