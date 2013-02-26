import orth
import seqs
import utils as ut
from pandas import DataFrame

def orth_count_ogroups(sp1, sp2):
    """
    symmetric measure of orthology
    """
    key, swap_order = orth.orth_key(sp1, sp2)
    ogs = orth._load_ogroups(ut.proj_path('convert_orth', 'table.'+key))
    return len(ogs)

def orth_count_genes(from_sp, to_sp, orth_func=orth.odict):
    """
    can use orth.odict_1to1 for symmetric measure
    """
    od = orth_func(from_sp, to_sp)
    return len(od)

def orth_count_table(sps, count_func=orth_count_genes, **kwargs):
    result = [[count_func(sp1, sp2, **kwargs) if sp2!=sp1 else count_genes(sp1) 
        for sp2 in sps] for sp1 in sps]
    return DataFrame(result, index=sps, columns=sps)

def count_genes(sp):
    ps = seqs.load_prots_from_fasta('../../data/sequences/canon/%s.fasta' % sp)
    return len(ps)


