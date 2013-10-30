from __future__ import division
import difflib
import itertools as it
import orth
import pairdict as pd
import seqs
import utils as ut

def load_paralogs(sp_base, sp_other):
    """
    Get all pairwise base species paralogs from sharing orthogroups in the
    inParanoid orthology files between the base and other species.
    """
    ogs = orth.load_ogroups(sp_base, sp_other)
    base_ogs = ut.i0(ogs)
    base_paralogs = groups_to_pairs(base_ogs)
    return base_paralogs

def groups_to_pairs(lol_groups):
    raw_pairs = ut.flatten([[x for x in it.combinations(group,2)] 
        for group in lol_groups])
    deduped = pd.dedupe(raw_pairs)
    return deduped

def difflib_ratio(aseq, bseq):
    """
    Given two protein sequences, compute the % identity between them.
    """
    diff = difflib.SequenceMatcher(a=aseq,b=bseq)
    return diff.ratio()

def identity_pairs(pairs, sp, compare_func=difflib_ratio):
    d_seqs = seqs.load_seqs_from_fasta(seqs.fasta_fname(sp))
    ratios = [(a, b, compare_func(d_seqs[a], d_seqs[b])) for a,b in pairs]
    return ratios

def paralog_identities(sp_base, sp_other):
    pairs = load_paralogs(sp_base, sp_other)
    idents = identity_pairs(pairs, sp_base)
    return idents

