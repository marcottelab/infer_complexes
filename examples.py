from __future__ import division
import itertools
from math import factorial
import numpy as np
import utils as ut
import random
from Struct import Struct
from pairdict import PairDict
import compare as cp

def base_examples_single(ppi_cxs, clean_cxs, all_cxs, split_fracs,
        both_cx_splits=None, ind_cycle=None):
    """
    Builds the training/test examples struct ready for scoring and learning.
    Test_neg_set: full list of proteins found in our data to use for test set
                  negatives generation. Use None to make negs the same way as
                  training, from the members of those complexes.
    ind_cycle: [0,-1] to start with train, then holdout, then test.  None for
    random sampling according to splits.
    """
    # Still not quite right: bashes any same-named complexes.
    # But 20121005 none of these exist currently.
    ppi_cxs,clean_cxs,all_cxs = [dict(cs) for cs in ppi_cxs,clean_cxs,all_cxs]
    lp_splits,clean_splits = (both_cx_splits or
            positives_from_corum(ppi_cxs, clean_cxs, split_fracs, ind_cycle))
    ptrain_lp = lp_splits[0]
    all_pos_lp = _complexes_to_LPairset(all_cxs)
    train_lp = add_negs(ptrain_lp, all_pos_lp, None, None)
    pdtrain = PairDict([(p[0],p[1],1 if p[2]=='true' else 0) 
        for p in train_lp.pairs]) 
    splits = (lp_splits, clean_splits)
    return pdtrain, splits

def base_examples(ppi_cxs, clean_cxs, all_cxs, splits, test_neg_set=None,
        nratio_train=None, nratio_test=None, pos_lengths=None, pos_splits=None,
        clean_splits=None, ind_cycle=None):
    """
    Builds the training/test examples struct ready for scoring and learning.
    Test_neg_set: full list of proteins found in our data to use for test set
                  negatives generation. Use None to make negs the same way as
                  training, from the members of those complexes.
    ind_cycle: [0,-1] to start with train, then holdout, then test.  None for
    random sampling according to splits.
    """
    # Still not quite right: bashes any same-named complexes.
    # But 20121005 none of these exist currently.
    ppi_cxs,clean_cxs,all_cxs = [dict(cs) for cs in ppi_cxs,clean_cxs,all_cxs]
    if pos_splits is None:
        pos_splits,clean_splits = positives_from_corum(ppi_cxs, clean_cxs,
              splits, ind_cycle)
    ptrain_lp,ptest_lp = pos_splits[:2]
    all_pos_lp = _complexes_to_LPairset(all_cxs)
    train_lp = add_negs(ptrain_lp, all_pos_lp, test_neg_set, nratio_train)
    test_lp = add_negs(ptest_lp, merge_lps([all_pos_lp,train_lp]),
            test_neg_set, nratio_test)
    return [PairDict([(p[0],p[1],1 if p[2]=='true' else 0) for p in lp.pairs])
            for lp in [train_lp, test_lp]], clean_splits


def lpairset_to_lol(lp):
    return list([list(tup) for tup in lp.pairs])

def add_negs(pos_lp, exclude_lp, from_set, ratio):
    """
    Add negatives in the ratio-to-postive specified by neg_ratios to the
    provided pos_splits, which is a list of LPairset--labeled pair set, ie a
    set of (ID1, ID2, 'true'/'false').  Take possible negatives for each split
    only from from_set, which for training is the positives set, and for
    testing is the set of all proteins with scores in our data.
    If a ratio isn't specified, just keep them all (asks for a huge number).
    """
    if from_set is None:
        from_set = pos_lp.members()
    k = len(pos_lp.pairs)*ratio if ratio is not None else 10**9
    nitems = len(from_set)
    if nitems*(nitems-1)/2 < 2*k:
        negs_lp = all_pairs(from_set, exclude_lp, k)
    else:
        negs_lp = kpairs(from_set, exclude_lp, k)
    full_lp = merge_lps([pos_lp, negs_lp])
    return full_lp

def all_pairs(items, exclude_lp, k):
    negs = [lpair(p1,p2,'false') for p1,p2 in itertools.combinations(
        items,2) if not exclude_lp.contains(lpair(p1, p2, 'true')) and not
        exclude_lp.contains(lpair(p1,p2,'false'))]
    return LPairset(set(maybe_sample(negs,k)))

def kpairs(items, exclude_lp, k):
    litems = list(items)
    n = 0
    selected = LPairset(set([]))
    while n<k:
        p1 = random.choice(litems)
        p2 = random.choice(litems)
        if (p1!=p2 and not exclude_lp.contains(lpair(p1,p2,'true')) and not
                exclude_lp.contains(lpair(p1,p2,'false'))):
            if selected.add(lpair(p1,p2,'false')):
                n = n+1
    return selected

def maybe_sample(pop,k,verbose=True):
    if k < len(pop):
        return random.sample(pop,k)
    else:
        if verbose:
            print "Can't sample: k", k, "larger than population", len(pop)
        return pop

def positives_from_corum(ppicxs, cleancxs, splits, ind_cycle,
                         split_for_unmatched=None):
    """
    Create training, [cross-val,] and test sets, and excluded set as the final
         one, with splits happening according to the supplied list of splits
                         and cleancxs complexes. Interactions forming the sets
                         come from ppicxs complexes.
    Ppicxs: dict{ complexA: set([p1,p2]); complexB: ...}
    Cleancxs: same format as ppicxs
    Splits: like [0,0.5,0.7,1] to form 50%,20%,30% approx subsets of
         interactions. [Training, Test, ..., Excluded]
    Assigns Ppicxs based on being subset of cleancxs; puts the few that are not
        subsets into splits[split_for_unmatched].
    """
    print "Total cleancxs:", len(cleancxs)
    if split_for_unmatched is None: 
        split_for_unmatched = random.choice(range(len(splits)-1))
    cleancxs2ppi = _merged_parents(ppicxs, cleancxs)
    clean_int_sets = [LPairset(set([lpair(x,y,'true') for p in cleancxs2ppi[c]
        for x,y in _set_to_pairs(ppicxs[p])]), name=(c, cleancxs[c]))
        for c in cleancxs2ppi]
    clean_ints = merge_lps(clean_int_sets)
    total_clean_pairs = len(clean_ints.pairs) # slow probably
    allsplits_lps, allsplits_cxs = positives_from_lps_sorted(clean_int_sets,
        total_clean_pairs, splits, ind_cycle)
    # Append the ppicxs interactions if any were left out.
    ppi_ints = _complexes_to_LPairset(ppicxs)
    diff_ppis_lp = ppi_ints.difference(clean_ints)
    if len(diff_ppis_lp.pairs)>0:
        assert False, "Ppi assignment problem: unmatched ppi cxs"
        #allsplits_lps[split_for_unmatched].merge(diff_ppis_lp)
        #allsplits_cxs[split_for_unmatched].append(str(len(diff_ppis_lp.pairs)) +
                                              #' ppi_cxs')
    allsplits_lps = dedupe_lps(allsplits_lps)
    print 'total clean ints:',total_clean_pairs
    print 'split interaction counts:', [len(x.pairs) for x in allsplits_lps]
    print 'split complex counts:', [len(x) for x in allsplits_cxs]
    print 'total assigned clean cxs:', sum([len(x) for x in allsplits_cxs])
    return allsplits_lps, allsplits_cxs

def _complexes_to_LPairset(cxs):
    lps_list = [LPairset(set([lpair(x,y,'true') 
        for x,y in _set_to_pairs(cxs[c])]),name=c) 
        for c in cxs]
    return merge_lps(lps_list)

def positives_from_lps_sorted(lps, total_pairs, splits, ind_cycle=None):
    """
    Random assignment to splits according to specified probabilities if
    ind_cycle is None; if ind_cycle specified, rotates according to that.
    """
    splits_lps = [LPairset(set([])) for s in splits[:-1]]
    splits_cxs = [[] for s in splits[:-1]]
    lps.sort(key=lambda lp: len(lp.pairs),reverse=True)
    index = ind_cycle[0] if ind_cycle else sample_index(splits)
    for lp in lps:
        split_lp = splits_lps[index]
        split_lp.merge(lp)
        splits_cxs[index].append(lp.name)
        index = ((index + ind_cycle[1]) % len(splits_lps) if ind_cycle else
                sample_index(splits))
    print "split positives count:", [len(slp.pairs) for slp in splits_lps]
    return splits_lps, splits_cxs

def sample_index(splits):
    r = random.random()
    for i,s in enumerate(splits[1:]):
        if r<s:
            return i

class LPairset(object):
    """
    Labeled Pair Set
    self.pairs: set of tuples of (id1,id2,'true/false')
    """

    def __init__(self, pairset, name=None):
        self.pairs = self.dedupe(pairset)
        if name: self.name = name

    def add(self, lpair):
        if not self.contains(lpair):
            self.pairs.add(lpair)
            return True # indicate success or not

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
        return LPairset(diff_set)

    def members(self):
        # returns a set of all members in all pairs
        return set.union(set([p[0] for p in self.pairs]),
                         set([p[1] for p in self.pairs]))

    def contains(self, lpair):
        return (lpair in self.pairs or lpair_flip(lpair) in self.pairs)

    def intersection(self, other):
        # returns pairs
        newlist = []
        for trip in self.pairs:
            if other.contains(trip):
                newlist.append(trip)
        return LPairset(set(newlist))
    
        
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
    # Changed 5/27: added first step to even out losses in training/test.
    newlps = [LPairset(lp.pairs) for lp in lps]
    def dedupe_fairly(lp1,lp2):
        intersect01 = lp1.intersection(lp2)
        int_half0 = LPairset(set(random.sample(list(intersect01.pairs),
                                               int(len(intersect01.pairs)/2))))
        int_half1 = intersect01.difference(int_half0)
        newlp1,newlp2 = lp1.difference(int_half0),lp2.difference(int_half1)
        return newlp1, newlp2
    for i,j in itertools.combinations(range(len(newlps)),2):
        #new1,new2 = dedupe_fairly(lps[i],lps[j]) WRONG 9/27. Re-adds dupes.
        newlps[i],newlps[j] = dedupe_fairly(newlps[i],newlps[j])
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
    print "len merged", len(merged)
    m2u = {}
    mcxs_names,mcxs_sets = zip(*merged.items())
    for uname,u in unmerged.items():
        mind = cp.best_match_index(u, mcxs_sets, cp.sensitivity)
        m2u.setdefault(mcxs_names[mind], set([])).add(uname)
    print "len m2u", len(m2u)
    return m2u

def _count_ints(nprots):
    def n_choose_r(n,r):
        return int(factorial(n) / ( factorial(r) * factorial(n-r) ))
    return n_choose_r(nprots,2) if nprots > 1 else 0
