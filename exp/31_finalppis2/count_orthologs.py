import elution as el
import orth
import seqs
import utils as ut
from pandas import DataFrame

def orth_count_ogroups(sp1, sp2):
    """
    Symmetric measure of orthology.
    Does not lend itself as well to only counting genes in a provided list.
    """
    key, swap_order = orth.orth_key(sp1, sp2)
    ogs = orth._load_ogroups(ut.proj_path('convert_orth', 'table.'+key))
    return len(ogs)

def orth_count_genes(from_sp, to_sp, orth_func=orth.odict, from_sp_genes=None):
    """
    orth.odict is asymmetric: any genes in from_sp with ortholog(s) in to_sp.
    Can use orth.odict_1to1 for symmetric measure.
    Checks whether the from_sp gene is in the provided set, from_sp_genes.
    If the same species given as both species, use all genes from that species.
    """
    if from_sp==to_sp:
        orths = all_genes(from_sp)
    else:
        orths = set(orth_func(from_sp, to_sp).keys())
    return len([g for g in orths if (from_sp_genes is None or g in from_sp_genes)])

def orth_count_table(sps, count_func=orth_count_genes, dict_sps_genes=None,
        **kwargs):
    result = [[count_func(sp1, sp2, from_sp_genes=dict_sps_genes[sp1] if
        dict_sps_genes else None, **kwargs) for sp2 in sps] for sp1 in sps]
    return result

def all_genes(sp):
    ps = seqs.load_prots_from_fasta('../../data/sequences/canon/%s.fasta' % sp)
    return ps


def df(result, sps):
    return DataFrame(data=result, index=sps, columns=sps)

def manysp_all_prots(sps, elutfs, **kwargs):
    d_allprots = dict([(s,el.all_prots(elutfs, sp_base=s, **kwargs)) 
        for s in sps])
    return d_allprots

