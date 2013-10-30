from __future__ import division
import elution as el
import orth
import seqs
import utils as ut
from pandas import DataFrame
from pandas import ExcelWriter

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


def df_sps(result, sps):
    return DataFrame(data=result, index=sps, columns=sps)

def set_from_pairs(ppis):
    return set(ut.i0(ppis) + ut.i1(ppis))

def manysp_all_prots(sps, elutfs, **kwargs):
    d_allprots = dict([(s,el.all_prots(elutfs, sp_base=s, **kwargs)) 
        for s in sps])
    return d_allprots

def found_orth_prots(base_set_prots, sp, base_sp):
    def isfound(d, key, set_values):
        for v in d.get(key, []):
            if v in set_values:
                return True
        return False
    od = orth.odict(sp, base_sp)
    if sp != base_sp:
        return set([g for g in all_genes(sp) if isfound(od, g, base_set_prots)])
    else:
        return base_set_prots

def manysp_protsets(sps, elutfs, lists_ppis, base_sp):
    """
    Orthology sharing computed from several sets of proteins:
    - Fasta files
    - In elution files, spectral count >=1
    - " " spectral count >=2
    - In each provided list of ppis [top_ppis, cxppis]
    """
    sp_prot_dicts = [
        dict([(sp, set(all_genes(sp))) for sp in sps]), 
        manysp_all_prots(sps, elutfs, min_count=1),
        manysp_all_prots(sps, elutfs, min_count=2)
        ]
    for ppis in lists_ppis:
        base_set_prots = set_from_pairs(ppis)
        dprots = dict([(sp, found_orth_prots(base_set_prots, sp, base_sp)) 
            for sp in sps])
        sp_prot_dicts.append(dprots)
    for d in sp_prot_dicts: 
        print "%s length:" % base_sp, len(d[base_sp])
    return sp_prot_dicts

def write_xls_sheets(dataframes, fname, sheet_names=None):
    writer = ExcelWriter(fname)
    sheet_names = sheet_names or ['sheet %s' %n for n in range(len(results))]
    for df,name in zip(dataframes, sheet_names):
        df.to_excel(writer, name)
    writer.save()


