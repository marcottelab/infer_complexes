from __future__ import division
import itertools
from math import factorial
import utils as ut
import random
from Struct import Struct

def base_examples(ppi_cxs, clean_cxs, splits=[0,0.33,0.66,1], neg_ratios=[4,40],
                  pos_lengths=None, shuf=True, pos_splits=None, confirm_redo=True):
    """
    Builds the training/test examples struct ready for scoring and learning.
    """
    if pos_splits is None:
        redo = True
        while redo:
            total_cleans = len(clean_cxs)
            pos_splits,ncleans = positives_from_corum(ppi_cxs, clean_cxs,
                  splits, shuf)
            #ensure we get a good number of complexes
            bal_fracs = [splits[i+1]-splits[i] for i in range(len(splits)-2)]
            balanced = min([len(x) > total_cleans*f for f,x in zip(bal_fracs,
                  ncleans[:len(neg_ratios)])]) 
            print 'balanced', [ total_cleans*f for f,x in zip(bal_fracs,
                  ncleans[:len(neg_ratios)])]
            redo = (shuf and not balanced) or (confirm_redo and
                  raw_input("Redo? y/n: ").lower() == 'y')
    full_splits = add_negs_to_splits(pos_splits, neg_ratios=neg_ratios)
    ltrain,ltest = [list([list(tup) for tup in s.pairs]) for s in full_splits]
    assert pos_lengths==None, 'Not implemented'
    for l in [ltrain,ltest]: random.shuffle(l)
    ex_struct = Struct(examples=ltrain+ltest, names=['id1','id2','hit'])
    return ex_struct, len(ltrain)

def add_negs_to_splits(pos_splits, neg_ratios):
    """
    Add negatives in the ratio-to-postive specified by neg_ratios to the
    provided pos_splits, which is a list of LPairset--labeled pair set, ie a
    set of (ID1, ID2, 'true'/'false').  Take possible negatives for each split
    only from the set of proteins with positive interactions in that split.
    """
    all_pos_lp = merge_lps(pos_splits)
    full_splits = []
    for lp,ratio in zip(pos_splits[:len(neg_ratios)],neg_ratios):
        members = lp.all_members()
        npos = len(lp.pairs)
        negs = []
        for i,p1 in enumerate(members):
            for j,p2 in enumerate(members):
                if j>i:
                    if not all_pos_lp.contains(lpair(p1,p2,'true')):
                        negs.append(lpair(p1,p2,'false'))
        sample_negs = len(lp.pairs)*ratio
        if sample_negs < len(negs):
            negs = random.sample(negs, sample_negs)
        else:
            print 'only', len(negs), 'negs, not', sample_negs
        full_splits.append(merge_lps([lp, LPairset(set(negs))]))
    return full_splits
    
def positives_from_corum(ppicxs, cleancxs, splits, shuffle_splits,
                         split_for_unmatched=1):
    """
    Create training, [cross-val,] and test sets, with splits happening
    according to the supplied list of splits and cleancxs complexes.
    Interactions forming the sets come from ppicxs complexes.
    Ppicxs: dict{ complexA: set([p1,p2]); complexB: ...}
    Cleancxs: same format as ppicxs
    Splits: like [0,0.5,0.7,1] to form 50%,20%,30% approx subsets of interactions.
    Assigns Ppicxs based on being subset of cleancxs; puts the few that are not
        subsets into splits[split_for_unmatched].
    """
    ppicxs = _remove_singles(ppicxs)
    cleancxs = _remove_singles(cleancxs)
    cleancxs2ppi = _merged_parents(ppicxs, cleancxs)
    clean_int_sets = [LPairset(set([lpair(x,y,'true') for p in cleancxs2ppi[c]
        for x,y in _set_to_pairs(ppicxs[p])]), name=c) for c in cleancxs2ppi]
    if shuffle_splits: random.shuffle(clean_int_sets)
    clean_ints = merge_lps(clean_int_sets)
    total_clean_pairs = len(clean_ints.pairs) # slow probably
    allsplits_lps, allsplits_cxs = positives_from_lps(clean_int_sets,
        total_clean_pairs, splits)
    # Append the ppicxs interactions
    ppi_int_sets = [LPairset(set([lpair(x,y,'true') for x,y in
        _set_to_pairs(ppicxs[p])]),name=p) for p in ppicxs]
    ppi_ints = merge_lps(ppi_int_sets)
    # Take the difference before merging, otherwise most ints end up here.
    allsplits_lps[split_for_unmatched].merge(LPairset(ppi_ints.difference(clean_ints)))
    allsplits_cxs[split_for_unmatched].append('ppi_cxs')
    allsplits_lps = dedupe_lps(allsplits_lps)
    print 'total clean ints:',total_clean_pairs
    print 'split interaction counts:', [len(x.pairs) for x in allsplits_lps]
    print 'split complex counts:', [len(x) for x in allsplits_cxs]
    return allsplits_lps, allsplits_cxs

def positives_from_lps(lps, total_pairs, splits):
    split_lp = LPairset(set([]))
    split_cxs = []
    allsplits_lps = []
    allsplits_cxs = []
    i = 0 # index for splits
    for lp in lps:
        threshold = total_pairs * (splits[i+1]-splits[i]) 
        maybe_split_lp = merge_lps([split_lp, lp])
        maybe_len = len(maybe_split_lp.pairs)
        if ( maybe_len > threshold and i < len(splits)-2 ):
            # if we've exceeded, only include this if it's mostly inside
            keep_it = (( maybe_len - threshold ) / len(lp.pairs) < 0.5 )
            if keep_it:
                allsplits_lps.append(maybe_split_lp)
                split_lp = LPairset(set([]))
                allsplits_cxs.append(split_cxs+[lp.name])
                split_cxs = []
            # otherwise save it for the next split
            else:
                allsplits_lps.append(split_lp)
                split_lp = lp
                allsplits_cxs.append(split_cxs)
                split_cxs = [lp.name]
            i = i+1
        else:
            split_lp = maybe_split_lp
            split_cxs.append(lp.name)
        #print i, 'ints', len(split_ints), [len(_dedup_lpair(x)) for x in \
        #ints_splits], 'cleancxs', len(split_cleancxs), [len(x) for x in cleancxs_splits]
    # Append the final split
    allsplits_lps.append(split_lp)
    allsplits_cxs.append(split_cxs)
    return allsplits_lps, allsplits_cxs

class LPairset(object):
    # Labeled Pair Set

    def __init__(self, pairset, name=None):
        self.pairs = self.dedupe(pairset)
        if name: self.name = name

    def add(self, lpair):
        if not self.contains(lpair):
            self.pairs.add(lpair)

    def merge(self, other):
        self.pairs = self.dedupe(set.union(self.pairs, other.pairs))

    def dedupe(self,pairset):
        # returns pairs
        new_set = pairset.copy()
        for trip in pairset:
            #if lpair_flip(trip) in new_set: much slower somehow
            if (trip[1],trip[0],trip[2]) in new_set:
                new_set.remove(trip)
        return new_set

    def difference(self, other):
        # returns pairs
        diff_set = self.pairs.copy()
        for trip in other.pairs:
            if trip in diff_set:
                diff_set.remove(trip)
            if lpair_flip(trip) in diff_set:
                diff_set.remove(lpair_flip(trip))
        return diff_set

    def all_members(self):
        # returns a list of all members in all pairs
        return set.union(set([p[0] for p in self.pairs]),
                         set([p[1] for p in self.pairs]))

    def contains(self, lpair):
        return (lpair in self.pairs or lpair_flip(lpair) in self.pairs)
    
        
    
        
def merge_lps(lps):
    def fmerge(first, second):
        newlps = LPairset(first.pairs)
        newlps.merge(second)
        return newlps
    return reduce(fmerge, lps)

def dedupe_lps(lps):
    # Cycles first through new, starting with new[1]
    # Must have only pairwise combinations instead of just i!=j
    # Otherwise you take away too many--remove overlapping pairs from both,
    # rather than from only one of the two new sets.
    newlps = [LPairset(lp.pairs) for lp in lps]
    for i,j in itertools.combinations(range(len(newlps)),2):
        # if random.random() > 0.5:
        #     x = i
        #     i = j
        #     j = x
        newlps[j] = LPairset(newlps[j].difference(lps[i]))
    # for i,_ in enumerate(lps):
    #     for j,_ in enumerate(newlps):
    #         if j>i:
    #             newlps[j] = LPairset(newlps[j].difference(lps[i]))
    return newlps

def lpair(x,y,z):
    return (x,y,z)

def lpair_flip(lpair):
    return (lpair[1],lpair[0],lpair[2])
    
def _set_to_pairs(s):
    return [(a,b) for a in s for b in s if a!=b]

    
    
def _remove_singles(complexes):
    return dict([(c,ps) for c,ps in complexes.items() if len(ps)>1])
    
def _merged_parents(unmerged, merged):
    """
    For now just return the first matching merged complex, based on subsets.
    Returns dict of merged_name:unmerged_name
    """
    def _complex_parents(u, unmerged, merged):
        return [k for k in merged.keys() if
            set.issubset(unmerged[u],merged[k])]
    m2u = []
    for u in unmerged:
        parents = _complex_parents(u,unmerged,merged)
        if len(parents)>0:
            m2u.append((parents[0], u))
    return ut.dict_sets_from_tuples(m2u)

def _count_ints(nprots):
    def n_choose_r(n,r):
        return int(factorial(n) / ( factorial(r) * factorial(n-r) ))
    return n_choose_r(nprots,2) if nprots > 1 else 0
