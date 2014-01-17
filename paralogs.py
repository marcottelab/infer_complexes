from __future__ import division
from collections import defaultdict
import itertools as it
import orth
import pairdict as pd
import ppi
import ppi_utils as pu
import seqs
import time
import utils as ut

def load_paralogs(sp_base, sp_other, ogroup_max, other_ogroup_max):
    """
    Get all pairwise base species paralogs from sharing orthogroups in the
    inParanoid orthology files between the base and other species.
    """
    ogs = orth.load_ogroups(sp_base, sp_other)
    use_ogs = []
    if ogroup_max:
        for og_base, og_other in ogs:
            if (len(og_base) <= ogroup_max 
                    and (other_ogroup_max is None
                        or len(og_other) <= other_ogroup_max)):
                use_ogs.append(og_base)
    else:
        use_ogs = ut.i0(base_ogs)
    base_paralogs = pu.groups_to_pairs(use_ogs)
    return base_paralogs

#def difflib_ratio(aseq, bseq):
    #"""
    #Given two protein sequences, compute the % identity between them.
    #This one seems like it really sucks.
    #"""
    #diff = difflib.SequenceMatcher(a=aseq,b=bseq)
    #return diff.ratio()

def identity_pairs(pairs, sp, compare_func=seqs.percent_identity, verbose=False):
    d_seqs = seqs.load_seqs_from_fasta(seqs.fasta_fname(sp))
    ratios=[]
    start = time.time()
    for i,(a,b) in enumerate(pairs):
        if verbose and (i % 10 == 0):
        #if verbose:
            mins = round((time.time() - start) / 60, 1)
            ut.printnow("%s of %s pairs, %s min" % (i, len(pairs), mins))
        try:
            ident = compare_func(d_seqs[a], d_seqs[b])
        except:
            print "too long or error, %s" % i 
            ident = -1
        ratios.append((a, b, ident))
    print "%s total not scored" % len([x for x in ratios if x[2]==-1])
    return ratios

def paralog_identities(sp_base, sp_other, ogroup_max=9, other_ogroup_max=None,
        verbose=False):
    pairs = load_paralogs(sp_base, sp_other, ogroup_max, other_ogroup_max)
    idents = identity_pairs(pairs, sp_base, verbose=verbose)
    return idents

def partner_scores_dict(ppis, setps):
    """
    Return { id1: [(id10, .5), (id23, .4)], id2: [(...), ...], ...]
    Process the full ppis list once to include all ids in the set setps.
    Will include asymmetric matches--must check when using for base-target
    interactions and exclude them.
    """
    assert type(setps) == set
    d = defaultdict(list)
    for ppi in ppis:
        if ppi[0] in setps:
            base = ppi[0]
            target = ppi[1]
        elif ppi[1] in setps:
            base = ppi[1]
            target = ppi[0]
        else:
            continue
        d[base].append((target, ppi[2]))
    return d

def interaction_distances(pairs, ppis):
    """
    For every pair in pairs, calculate the distance between their interaction
    patterns.
    """
    allps = pu.unique_items(pairs)
    dpartners = partner_scores_dict(ppis, allps)
    dists = []
    for a,b in pairs:
        dists.append(distance_partners(a, b, dpartners))


def distance_partners(a, b, dpartners):
    #aints,bints = [get_partner_scores(x, ppis) for x in a,b]
    dcomb = combine_scores(a, b, dpartners)
    return dcomb

def combine_scores(a, b, dpartners, default=0.01):
    """
    Match up the partners eg { id1: [ascore, bscore], id2: [default, bscore]}
    """
    d = defaultdict(lambda: [default, default])
    for pbase, index in [(a, 0), (b, 1)]:
        for target,score in dpartners[pbase]:
            d[target][index] = score
    return d

#def combine_scores_ppis(a, b, ppis, default=0.01):
    #"""
    #Match up the partners eg { id1: [ascore, bscore], id2: [default, bscore]}
    #"""
    #d = defaultdict(lambda: [default, default])
    #for pbase,index in [(a, 0), (b, 1)]:
        #for ppi in ppis:
            #if ppi[0]==pbase:
                #target = ppi[1]
            #elif ppi[1]==pbase:
                #target = ppi[0]
            #else:
                #continue
            #d[target][index] = ppi[2]
    #return d

def get_partner_scores(base, ppis):
    accum = []
    for p in ppis:
        if p[0]==base:
            accum.append((p[1],p[2]))
        elif p[1]==base:
            accum.append((p[0],p[2]))
    return accum

def score_pairs_nonpairs(pairs, elut_fs, spbase, scores=['poisson'],
        score_cutoff=0.25, do_nonpairs=False):
    pdpairs = pd.PairDict(pairs)
    for p in pdpairs.d: pdpairs.d[p] = 1 #marked as true paralog
    if do_nonpairs:
        for np in pu.nonpairs_gen(pairs, len(pairs)):
                pdpairs.set(np, 0) # marked as non-paralog
    arr = ppi.new_score_array(pdpairs, scores, elut_fs, [])
    del pdpairs #lots of memory
    scored_arr = ppi.score_and_filter(arr, scores, elut_fs, score_cutoff,
            spbase, [], '', do_filter=False)
    return scored_arr
  
