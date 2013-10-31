from __future__ import division
import difflib
import itertools as it
import orth
import pairdict as pd
import ppi_utils as pu
import seqs
import utils as ut

def load_paralogs(sp_base, sp_other, ogroup_limit):
    """
    Get all pairwise base species paralogs from sharing orthogroups in the
    inParanoid orthology files between the base and other species.
    """
    ogs = orth.load_ogroups(sp_base, sp_other)
    base_ogs = ut.i0(ogs)
    if ogroup_limit:
        base_ogs = [x for x in base_ogs if len(x) < ogroup_limit]
    base_paralogs = pu.groups_to_pairs(base_ogs)
    return base_paralogs

def difflib_ratio(aseq, bseq):
    """
    Given two protein sequences, compute the % identity between them.
    This one seems like it really sucks.
    """
    diff = difflib.SequenceMatcher(a=aseq,b=bseq)
    return diff.ratio()

def identity_pairs(pairs, sp, compare_func=seqs.percent_identity):
    d_seqs = seqs.load_seqs_from_fasta(seqs.fasta_fname(sp))
    ratios = [(a, b, compare_func(d_seqs[a], d_seqs[b])) for a,b in pairs]
    return ratios

def paralog_identities(sp_base, sp_other, ogroup_limit=10):
    pairs = load_paralogs(sp_base, sp_other, ogroup_limit)
    idents = identity_pairs(pairs, sp_base)
    return idents

