from __future__ import division
from collections import defaultdict
import corum as co
import numpy as np
import pairdict as pd
import ppi_utils as pu
import score as sc
from Struct import Struct
import utils as ut

def load_ppis(fname='/Users/blakeweb/Dropbox/complex/shared_complexes/predictions/Hs_2sp_35.6_ppis_highconfidence_from35.3_ppis_goldlabeled_top16655_cumprec70pct.txt'):
    return pu.load_ppis(fname)

def load_cxppis(fname='/Users/blakeweb/Dropbox/complex/shared_complexes/predictions/Hs_2sp_v35.6_clustered_ppis_981cxs.txt'):
    return pu.load_ppis(fname)

def load_hpa_localization(fname='./enrichment_datasets/subcellular_location.csv'):
    locs = ut.load_lot(fname, sep=",")[1:]
    locs_clean = [[x.strip("\"") for x in l] for l in locs]
    locs_filt = [l for l in locs_clean if l[4]=="Supportive"]
    locs = [(l[0],l[1]) for l in locs_filt]
    return locs

def attr_to_sets(atts):
    """
    Given a list of items with attributes, return a list of item sets.
    Input: [(item1, attr1), (item2, attr2), (item3, attr1), ...]
    Output: [{item1,item3}, {item2}, ...]
    """
    datts = defaultdict(set)
    for i,a in atts:
        datts[a].add(i)
    return ut.i1(datts.items())

def items_from_pairs(pairs):
    return set(ut.i0(pairs)+ut.i1(pairs))

def pair_intersect(a, b):
    pda = pd.PairDict(a)
    return len([p for p in b if pda.contains(p)])

def filter_pairs(pairs, items):
    items = set(items)
    return [p for p in pairs if p[0] in items and p[1] in items]

def hpa_stats(ppis, locs, max_set_size=None):
    s = attr_to_sets(locs)
    if max_set_size is not None: 
        s = [c for c in s if len(c) < max_set_size]
    plocs = co.pairs_from_complexes(s)
    ppiprots = set(ut.i0(ppis)+ut.i1(ppis))
    anprots = set(ut.i0(locs))
    intprots = set.intersection(ppiprots, anprots)
    print len(ppiprots), len(anprots), len(intprots)
    return ppis_stats(ppis, plocs, intprots)

def ppis_stats(ppis, ppis_benchmark, bg_prots):
    """
    all_prots_1: set of observable proteins forming our predicted ppis, 
    """
    #bg_prots = set.intersection(items_from_pairs(ppis_benchmark), all_prots)
    print len(bg_prots)
    n_bg_ppis = int((len(bg_prots) * (len(bg_prots)-1)) / 2)
    filt_ppis = filter_pairs(ppis, bg_prots)
    filt_bench = filter_pairs(ppis_benchmark, bg_prots)
    n_intersect = pair_intersect(filt_ppis, filt_bench)
    k, N, m, n = (n_intersect, n_bg_ppis, len(filt_ppis), len(filt_bench))
    print "(for hyperg_cum:)"
    print "k (intersection): %s; N (bg): %s; m (#ppis): %s; n: %s" %(k,N,m,n)
    n_filt_bench = len(filt_bench)
    print "Expected: %s of %s (%s%%)" % (n_filt_bench, n_bg_ppis,
            np.round(100*n_filt_bench/n_bg_ppis, 2))
    n_annprots_ppis = len(filt_ppis)
    print "Observed: %s of %s (%s%%)" % (n_intersect, n_annprots_ppis,
            np.round(100*n_intersect/n_annprots_ppis, 2))
    #Observed: %s; Fold change: %s
    #print "Min set size:", min(len(ppis1), len(ppis2))
    #hyper = hyperg_cum(intersect, total_ppis, len(ppis1), len(ppis2),
            #verbose=True)
    #return intersect, hyper
    return k,N,m,n
