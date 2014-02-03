from __future__ import division
import utils as ut
from collections import defaultdict
from collections import Counter

def sums_proteins_ppis(ppis):
    """
    Return a ranked list of [(geneid, sum of ppis) .. ]
    """
    ppis = [(p[0],p[1],float(p[2]),int(p[3])) for p in ppis]
    psis = [(p[0],p[2]) for p in ppis] + [(p[1],p[2]) for p in ppis]
    dsums = defaultdict(float)
    for p in psis:
        dsums[p[0]] += p[1]
    sums = dsums.items()
    sums.sort(key=lambda x: x[1], reverse=True)
    return sums

def sums_proteins_complexes(cxs):
   counts = Counter(ut.flatten(cxs)).items()
   counts.sort(key=lambda x: x[1], reverse=True)
   return counts
    


