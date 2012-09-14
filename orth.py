import utils as ut
import os

keys = "Hs-Ce Hs-Dd Hs-Dm Hs-Mm Hs-Sp Hs_uni-Ce_uni Ce-Dm Ce-Mm Ce-Sp Sp-Dm Sp-Mm".split()

def odict(from_sp, to_sp):
    """
    Load a dict from file, eg:
    {HsProt1: set([CeProtA, CeProtB,...]), ...}
    """
    key = from_sp + '-' + to_sp
    if key in keys:
        swap_order=False
    else:
        key = to_sp + '-' + from_sp
        if key in keys:
            swap_order=True
        else:
            assert False, "Orthogroup key not in keys list"
    return _ogroups_to_odict(_load_ogroups(ut.proj_path('convert_orth',
                                                        'table.'+key)),
                             swap_order=swap_order)

def convert_dict_single(fromtype, totype):
    """
    totype: must be Sp (eg 'Hs') or Sp_seqdb
    Returns None if not necessary or not found.
    """
    if len(totype.split('_')) > 1:
        # Get rid of the 2nd half of totype if it's default for that species
        tosp, toseqdb = totype.split('_')
        if toseqdb == ut.config()[tosp+'_default']:
            totype = tosp
    if fromtype == totype:
        return None
    elif len(fromtype) == len(totype) == 2:
        return odict(fromtype, totype)
    else:
        return custom_conversion(fromtype, totype)

def convert_dict(fromtype, totype):
    """
    First looks for single conversion step. If not found, splits it up.
    Returns None if not necessary or not found.
    """
    conv1 = convert_dict_single(fromtype, totype)
    if conv1:
        return conv1
    else:
        # If we made it here, try first converting to second species, 
        # then looking for other conversion.
        conv1 = convert_dict_single(fromtype, totype[:2])
        conv2 = convert_dict_single(totype[:2], totype)
        if conv1 and conv2:
            return ut.compose_dict_sets(conv1,conv2)

def custom_conversion(fromtype, totype):
    """
    Check for a custom file in data/convert
    Return None if not found.
    """
    fname = "%s2%s.tab" % (fromtype, totype)
    fpath = ut.proj_path('convert',fname)
    if os.path.exists(fpath):
        return ut.load_dict_sets(fpath)


            
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
