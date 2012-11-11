from __future__ import division
import itertools as it
from collections import defaultdict
import numpy as np
import networkx as nx
import utils as ut
import pairdict as pd

def poss_ints(nps):
    return (nps*(nps-1))/2

def cliqueness_mult(cxs,cxppis):
    pdppis = pd.PairDict(cxppis)
    return [cliqueness(c, pdppis) for c in cxs]

def cliqueness(cx,pdppis):
    return len([1 for p1,p2 in it.combinations(cx,2) if
        pdppis.contains((p1,p2))]) / poss_ints(len(cx)) 

def filter_cxs(cxs, cxppis, qness=.7, min_size=0):
    pdppis = pd.PairDict(cxppis)
    return [c for c in cxs if cliqueness(c,pdppis)>qness and len(c)>=min_size]

def nx_cliques(ppis, min_len=3, min_weight=0):
    G = nx.Graph()
    G.add_weighted_edges_from([p[:3] for p in ppis])
    qs = [set(c) for c in nx.find_cliques(G) if len(c) >= min_len]
    if min_weight:
        qs = [q for q in qs if avg_weight(G,q) > min_weight]
    return qs

def avg_weight(G, clique):
    return np.mean([G[p1][p2]['weight'] 
        for p1,p2 in it.combinations(clique,2)])

def genes2phens(p2entrez_file, e2z=None):
    ensg2z = e2z if e2z else ut.load_dict_sets(ut.proj_path('convert', 'Hs2Hs_entrez.tab'))
    def dict_sets_rem_prefix(d, sep):
        d = dict([(k,set([vi.split(sep)[1] for vi in v])) for k,v in d.items()])
        return d
    p2z = dict_sets_rem_prefix(ut.load_dict_sets(p2entrez_file), ":")
    return ut.compose_dict_sets(ensg2z, ut.dict_inverse_sets(p2z))

