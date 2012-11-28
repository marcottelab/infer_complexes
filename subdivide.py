from __future__ import division
import itertools as it
import numpy as np
import utils as ut
import evidence as ev
import ml
import pairdict as pd
import elution as el

empty_val = 'Empty'

class SubDivider(object):

    def __init__(self, arr, cluster, elutfs, pd_feats=None, elmax=None,
            avg_ev=0.2):
        self.pd_clust_feats, self.feat_names = (pd_evidences(cluster, arr) if pd_feats is None
                else pd_feats)
        self.elut_max = (elut_gene_maxes(elutfs, cluster) if elmax is None else
                elmax)
        self.cluster = cluster
        self.avg_ev = avg_ev

    def go_split(self):
        # all combs of 1/0 assignments for each member of the cluster
        n = len(self.cluster)
        all_ins_outs = [np.binary_repr(x,width=n) for x in range(2**n)]
        raw_scores = [self.score_split(zip(self.cluster, ins_outs)) for
                ins_outs in all_ins_outs]
        print raw_scores
        scores = to_z_avg_wempties(raw_scores)
        scored_splits = zip(scores, all_ins_outs)
        return scored_splits

    def score_split(self, members):
        """
        members: list of [(id1, incomplex), ...]
        """
        return [self.score_pair(i,j) for i,j in it.combinations(members,2)]

    def score_pair(self, i_i_in, j_j_in):
        i,i_in = i_i_in; j,j_in = j_j_in
        i_in, j_in = int(i_in), int(j_in)
        if i_in and j_in:
            return self.score_together(i,j)
        elif i_in or j_in:
            return self.score_apart(i,j)
        else:
            return empty_val

    def score_together(self, i, j):
        """
        Looking only for evidence they're together in ANY contexts, regardles of
        wehther they're also not in others >> clip to include only positive values.
        """
        pair = self.pd_clust_feats.find((i,j))
        return sum(np.clip(self.pd_clust_feats.d[pair],0,1000)) if pair else 0


    def score_apart(self, gi, gj, cutoff=2):
        pair = self.pd_clust_feats.find((gi,gj))
        if pair: 
            evs = self.pd_clust_feats.d[pair]
        else:
            evs = [0]*len(self.feat_names)
        cumsum = 0
        for elutf in self.elut_max:
            maxi,maxj = [self.elut_max[elutf].get(x,0) for x in gi,gj]
            if maxi >= cutoff and maxj >= cutoff:
                use_inds = [ind for ind,name in enumerate(self.feat_names) if
                        frac_name(name) == ut.shortname(elutf)]
                cumsum += sum([self.avg_ev-evs[ind] for ind in use_inds])
            #elif maxi > cutoff or maxj > cutoff:
                #pass
        return cumsum
        #else:
            #return 0

def frac_name(score_name):
    return '_'.join(score_name.split('_')[:-1])

def to_z_avg_wempties(lol):
    def flat(lol):
        return [x for l in lol for x in l if x!=empty_val]
    mu,sd = np.mean(flat(lol)), np.std(flat(lol))
    def zscore_wempties(x, mu, sd):
        return 0 if x==empty_val else (x-mu)/sd
    newlol = [[zscore_wempties(x, mu, sd) for x in l] for l in lol]
    return [np.mean(l) for l in newlol]

def elut_gene_maxes(elutfs, geneids):
    d = {}
    for f in elutfs:
        e = el.load_elution(f)
        prots_inv = ut.list_inv_to_dict(e.prots)
        for gid in geneids:
            if gid in prots_inv:
                d.setdefault(f,{})[gid] = np.max(e.mat[prots_inv[gid]])
    return d

def pd_evidences(cluster,arr):
    arr_clust = arr[[i for i,r in enumerate(arr) 
            if r[0] in cluster and r[1] in cluster]]
    # doesn't seem right--if most the interactions are strong, shouldn't
    # normalize them down--should still count. but does it matter since it's
    # all comparative? messes with thinking about the clipping though in
    # score_together.
    #features, _ = ml.normalize(ml.arr_feats(arr_clust))
    sps = ut.config()['elut_species'].split('_')
    names = [n for n in arr.dtype.names[3:] if n[:2] in sps]
    features = arr[names]
    pd_ev = pd.PairDict([[r[0],r[1]] + list(features[i]) for i,r in
        enumerate(arr_clust)])
    return pd_ev, names
