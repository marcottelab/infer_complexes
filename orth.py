import utils as ut
import os

def odict(from_sp, to_sp, whichseqset='ensembl/longest'):
    """
    Load a dict from file, eg:
    {HsProt1: set([CeProtA, CeProtB,...]), ...}
    """
    key = from_sp + '-' + to_sp
    if key in keys:
        return _ogroups_to_odict(_load_ogroups(_inpar_file(key, whichseqset)))
    else:
        key = to_sp + '-' + from_sp
        if key in keys:
            return _ogroups_to_odict(_load_ogroups(_inpar_file(key, whichseqset)),
                    swap_order=True)
        else:
            assert False, "Orthogroup key not found"
            
def _inpar_file(key, whichseqset):
    return os.path.expanduser('~/Dropbox/complex/data/orthologies/' +
                whichseqset + '/table.' + key)
keys = "Hs-Ce Hs-Dd Hs-Dm Hs-Mm Hs-Sp Hs.fa-Ce.fa".split()

def _ogroups_to_odict(ogroups, swap_order=False):
    """
    From a list of orthogroups, return a dict from sp1 prots to a set of sp2
    prots. We want a dictionary from the first species in the file to the second,
    unless swap_order is True.  
    """
    sp1col = 1 if swap_order else 0
    sp2col = 0 if swap_order else 1
    orthdict = dict([(p1,set([p2 for p2 in og[sp2col]])) for og in ogroups for
                p1 in og[sp1col]])
    return orthdict

def _load_ogroups(fname):
    """
    Load an inparanoid table.Sp1-Sp2 file into a list of orthogroups, where
    each orthogroup is a tuple containing 1) a list of proteins in sp1 and 2) a
    list of proteins in sp2.
    Eg: [([HsProtA, HsProtB,..],[CeProtA,CeProtB,..]), ([..],[..]), ...]
    """
    # Skip header row; protein ids alternate with meaningless conf scores in
    # columns 2 and 3 in the order of the filename
    ogroups = [([p for p in row[2].split()[::2]],[p for p in
            row[3].split()[::2]]) for row in ut.load_tab_file(fname)][1:]
    return ogroups
