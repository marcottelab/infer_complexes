from __future__ import division
import itertools
from math import factorial
import utils as ut
import random

def add_negs_to_splits(pos_splits, neg_ratios=[1,10]):
    """
    Add negatives in the ratio-to-postive specified by neg_ratios to the
    provided pos_splits, which is a list of LPairset--labeled pair set, ie a
    set of (ID1, ID2, 'true'/'false').  Take possible negatives for each split
    only from the set of proteins with positive interactions in that split.
    """
    all_pos_lp = merge_lps(pos_splits)
    full_splits = []
    for lp,ratio in ut.zip_exact(pos_splits[:len(neg_ratios)],neg_ratios):
        members = lp.all_members()
        npos = len(lp.pairs)
        negs = []
        for i,p1 in enumerate(members):
            for j,p2 in enumerate(members):
                if j>i:
                    if not all_pos_lp.contains(lpair(p1,p2,'true')):
                        negs.append(lpair(p1,p2,'false'))
        negs = random.sample(negs, len(lp.pairs)*ratio)
        full_splits.append(merge_lps([lp, LPairset(set(negs))]))
    return full_splits
    
def positives_from_corum(ppicxs, cleancxs, splits=[0,0.3,0.6,1],
        cleancxs_shuffle=False, split_for_unmatched=1):
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
    clean_int_sets = [LPairset(set([lpair(x,y,'true') for p in cleancxs2ppi[c] for
        x,y in _set_to_pairs(ppicxs[p])]), name=c) for c in cleancxs]
    if cleancxs_shuffle: random.shuffle(clean_int_sets)
    clean_ints = merge_lps(clean_int_sets)
    total_clean_pairs = len(clean_ints.pairs) # slow probably
    print 'total clean ints:',total_clean_pairs
    allsplits_lps, allsplits_cxs = positives_from_lps(clean_int_sets, total_clean_pairs, splits)
    # Append the ppicxs interactions
    ppi_int_sets = [LPairset(set([lpair(x,y,'true') for x,y in
        _set_to_pairs(ppicxs[p])]),name=p) for p in ppicxs]
    ppi_ints = merge_lps(ppi_int_sets)
    # Take the difference before merging, otherwise most ints end up here.
    allsplits_lps[split_for_unmatched].merge(LPairset(ppi_ints.difference(clean_ints)))
    allsplits_cxs[split_for_unmatched].append('ppi_cxs')
    return dedupe_lps(allsplits_lps), allsplits_cxs

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
            if lpair_flip(trip) in new_set:
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

def pos_neg_pairs(complexes, complexes_exclude):
    """
    Complexes: a dict of all the complexes to use in this part of the process.
    If ppi, this should be the ppi-specific unmerged complexes.  If for complex
    prediction, this should be the complex-specific merged complexes.
    complexes_exclude: a dict of all the complexes whose interactions should be
    excluded from the learning set.
    This complexity is necessary to be able to
    use the different complex sets for ppi learning vs complex learning.
    """
    prots = list(reduce(set.union,[complexes[k] for k in complexes]))
    true_ints = corum_ints_duped(complexes)
    exclude_ints = corum_ints_duped(complexes_exclude)
    pos = []
    pos_ex = []
    negs = []
    for ind, i in enumerate(prots):
        for j in prots[ind:]:
            if i in true_ints and j in true_ints[i]:
                if i in exclude_ints and j in exclude_ints[i]:
                    pos_ex.append((i,j,'true'))
                else:
                    pos.append((i,j,'true'))
            else:
                if i not in exclude_ints or j not in exclude_ints[i]:
                    negs.append((i,j,'false'))
    for l in [pos,pos_ex,negs]: random.shuffle(l)
    plen = len(pos)
    split = int(len(negs) * (plen / (plen + len(pos_ex))))
    negs_use = negs[:split]
    negs_ex = negs[split:]
    for l,f in zip([pos,pos_ex,negs_use,negs_ex],[fname_pos,
            ut.pre_ext(fname_pos, '_exclude'), ut.pre_ext(fname_pos, '_negs'),
            ut.pre_ext(fname_pos, '_negs_exclude')]):
        ut.write_tab_file(l,f)

def full_examples(key, elut_fs, scores, species, fnet_gene_dict,
                  elut_score_cutoff=0.5):
    """
    Key like 'Ce_ensp', 'Hs_uni'. species like 'Hs'.
    Use fnet_gene_dict = -1 to skip functional network.  None means no dict is
        needed. Can supply the dict itself or a string--like 'cep2ceg' or
        'paper_uni2ensg'
    npos=nnegs=None means to use all pos and matching length negs.
    For the set of train_frac, load equal pos and neg.  For remainder (test)
        load nnegs negs.
    """
    # Train and test are merged then split to speed this up 2x
    ex_struct, ntrain = base_examples(key)
    el.score_multi_elfs(ex_struct, elut_fs, scores)
    # Filter out train AND test examples without a score exceeding cutoff
    if elut_fs and elut_score_cutoff is not None:
        ex_struct, ntrain = split_filt_merge(ex_struct, range(3,
                  len(ex_struct.names)), elut_score_cutoff, ntrain)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_struct, species, genedict=fnet_gene_dict)
    exs_train, exs_test = exstruct_split(ex_struct, ntrain)
    return exs_train, exs_test

def predict_all(elut_fs, scores, species, fnet_gene_dict, elut_score_cutoff=0.5):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    pairs = el.all_filtered_pairs(elut_fs, scores, elut_score_cutoff)
    # examples like [['id1', 'id2', 'true/false'], ...]
    exs = [[p1, p2, '?'] for p1,p2 in pairs]
    ex_struct = Struct(examples=exs,names=['id1','id2','hit'])
    el.score_multi_elfs(ex_struct, elut_fs, scores)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_struct, species, genedict=fnet_gene_dict)
    return ex_struct
