from __future__ import division
from math import factorial
import utils as ut

def fixit(ppicxs, cleancxs, splits=[0,0.3,0.6,1],
        cleancxs_shuffle=True, split_for_unmatched=1):
    """
    Create training, [cross-val,] and test sets, with splits happening
    according to the supplied list of splits and cleancxs complexes.
    Interactions forming the sets come from ppicxs complexes.
    Ppicxs: dict{ complexA: set([p1,p2]); complexB: ...}
    Cleancxs: same format as ppicxs
    Splits: like [0,0.5,0.7,1] to form 50%,20%,30% approx subsets of interactions.
    Assigns Ppicxs based on being subset of cleancxs; handles non subsets
    proportionally.
    """
    ppicxs = _remove_singles(ppicxs)
    cleancxs = _remove_singles(cleancxs)
    cleancxs2ppi = _merged_parents(ppicxs, cleancxs)
    ints_splits = [] # container for lists of interaction pairs
    cleancxs_splits = []
    clean_int_sets = [LPairset(set([lpair(x,y,'true') for p in cleancxs2ppi[c] for
        x,y in _set_to_pairs(ppicxs[p])]), name=c) for c in cleancxs]
    clean_ints = merge_lps(clean_int_sets)
    # print 'total ppicxs:', len(ppi_ints.pairs)
    # print 'total', len(merge_LPSets([clean_ints,ppi_ints]).pairs)
    # diff_ints = ppi_ints.difference(clean_ints)
    # print 'difference', len(diff_ints)
    if cleancxs_shuffle: random.shuffle(clean_int_sets)
    allsplits_lps, allsplits_cxs = positives_from_lps(clean_int_sets, splits)
    # Append the ppicxs interactions
    ppi_int_sets = [LPairset(set([lpair(x,y,'true') for x,y in
        _set_to_pairs(ppicxs[p])]),name=p) for p in ppicxs]
    ppi_ints = merge_lps(ppi_int_sets)
    allsplits_lps[split_for_unmatched].merge(ppi_ints)
    allsplits_cxs[split_for_unmatched].append('ppi_cxs')
    allsplits_lps = dedupe_lps(allsplits_lps)
    return allsplits_lps, allsplits_cxs

def positives_from_lps(lps, splits):
    total_pairs = len(merge_lps(lps).pairs) # slow probably
    print 'total clean ints:',total_pairs
    split_lp = LPairset(set([]))
    split_cxs = []
    allsplits_lps = []
    allsplits_cxs = []
    i = 0 # index for splits
    for lp in lps:
        threshold = total_pairs * (splits[i+1]-splits[i]) 
        maybe_split_lp = merge_lps([split_lp, lp])
        if ( len(maybe_split_lp.pairs) > threshold and i < len(splits)-2 ):
            # half the time include this last one
            allsplits_lps.append(maybe_split_lp)
            split_lp = LPairset(set([]))
            allsplits_cxs.append(split_cxs+[lp.name])
            split_cxs = []
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
        
def merge_lps(lps):
    def fmerge(first, second):
        newlps = LPairset(first.pairs)
        newlps.merge(second)
        return newlps
    return reduce(fmerge, lps)

def dedupe_lps(lps):
    # Cycles first through new, starting with new[1]
    # Must have the j>i constraint or its reverse instead of just i!=j
    # Otherwise you take away too many--remove overlapping pairs from both,
    # rather than from only one of the two new sets.
    newlps = [LPairset(lp.pairs) for lp in lps]
    for i,_ in enumerate(lps):
        for j,_ in enumerate(newlps):
            if j>i:
                newlps[j] = LPairset(newlps[j].difference(lps[i]))
    return newlps

def _flatten_ints(complex_ints):
    return _dedup_lpair(reduce(set.union,(ints for c,ints in
        complex_ints)))

def lpair(x,y,z):
    return (x,y,z)

def lpair_flip(lpair):
    return (lpair[1],lpair[0],lpair[2])
    
def _set_to_pairs(s):
    return [(a,b) for a in s for b in s if a!=b]

def _deep_dedup_lpair(trip_sets):
    # BUG: This is WRONG.
    # Taking away too many.  I end up removing any overlaps from both, rather
    # than one, of the new sets.
    new_sets = [s.copy() for s in trip_sets]
    for i,s in enumerate(trip_sets):
        for j,t in enumerate(new_sets):
            for trip in s:
                trip_rev = lpair(trip[1],trip[0],trip[2])
                if trip_rev in t:
                    t.remove(trip_rev)
                if i!=j and trip in t:
                    t.remove(trip)
    return new_sets
    
def _dedup_lpair(trips_set):
    new_set = trips_set.copy()
    for trip in trips_set:
        if lpair(trip[1],trip[0],trip[2]) in new_set:
            new_set.remove(trip)
    return new_set
    
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

def examples_from_complexes(unmerged, merged, splits=[0,0.3,0.6,1],
        merged_shuffle=True, split_for_unmatched=1):
    """
    Create training, [cross-val,] and test sets, with splits happening
    according to the supplied list of splits and merged complexes.
    Interactions forming the sets come from unmerged complexes.
    Unmerged: dict{ complexA: set([p1,p2]); complexB: ...}
    Merged: same format as unmerged
    Splits: like [0,0.5,0.7,1] to form 50%,20%,30% approx subsets of interactions.
    Assigns Unmerged based on being subset of merged; handles non subsets
    proportionally.
    """
    unmerged = _remove_singles(unmerged)
    merged = _remove_singles(merged)
    merged2un = _merged_parents(unmerged, merged)
    ints_splits = [] # container for lists of interaction pairs
    merged_splits = []
    merged_ints = [(m, set([lpair(x,y,'true') for u in merged2un[m]
        for x,y in _set_to_pairs(unmerged[u])])) for m in merged]
    unmerged_ints = [(u, set([lpair(x,y,'true') for u in unmerged
        for x,y in _set_to_pairs(unmerged[u])]))]
    total_merged = len(_dedup_lpair(_flatten_ints(merged_ints)))
    print 'total merged:',total_merged
    if merged_shuffle: random.shuffle(merged_ints)
    split_ints = set([])
    split_merged = []
    i = 0 # index for splits
    for m,new_ints in merged_ints:
        new_ints = _dedup_lpair(new_ints)
        new_split_ints = _dedup_lpair(set.union(split_ints, new_ints))
        if ( len(new_split_ints) > total_merged * (splits[i+1]-splits[i]) and
             i < len(splits)-2 ):
            if random.random() > 0.5:
                # half the time include this last one
                ints_splits.append(new_split_ints)
                split_ints = set([])
                merged_splits.append(split_merged+[m])
                split_merged = []
            else:
                ints_splits.append(split_ints)
                split_ints = new_ints
                merged_splits.append(split_merged)
                split_merged = [m]
            i = i+1
        else:
            split_ints = new_split_ints
            split_merged.append(m)
        print i, 'ints', len(split_ints), [len(_dedup_lpair(x)) for x in \
        ints_splits], 'merged', len(split_merged), [len(x) for x in merged_splits]
    # Append the final split
    ints_splits.append(split_ints)
    merged_splits.append(split_merged)
    # Append the unmerged interactions
    unmerged_to_add = _dedup_lpair(set.difference(_flatten_ints(unmerged_ints),
        _flatten_ints(merged_ints)))
    ints_splits[split_for_unmatched] = set.union(
        ints_splits[split_for_unmatched], unmerged_to_add)
    merged_splits[split_for_unmatched].append('Unmerged: %s interactions' %
        len(unmerged_to_add))
    return _deep_dedup_lpair(ints_splits), merged_splits
