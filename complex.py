import os
import utils as ut

corum_files = {
    'uppi': '~/Dropbox/complex/data/corum/corum_overlaps_forppi.tab',
    'eppi': '~/Dropbox/complex/data/corum/corum_overlaps_forppi_en.tab',
    'uclean': '~/Dropbox/complex/data/corum/corum_clean_nooverlaps_supptable3.tab',
    'eclean':
        '~/Dropbox/complex/data/corum/corum_clean_nooverlaps_supptable3_ensp.tab'
    }

def corum_ints_duped(key):
    return _load_intdict(os.path.expanduser(corum_files[key]))

def corum_pairs(key):
    return _load_complex_pairs(os.path.expanduser(corum_files[key]))

def _load_complex_pairs(filename):
    intdict = _load_intdict(filename)
    int_dedup = ut.dict_dedup(intdict)
    pairs = []
    for p, partners in int_dedup.items():
        for par in partners:
            pairs.append((p,par))
    return pairs

def load_complexes(filename, format_single_protein=False):
    # load corum-type file into a dictionary
    # complexes: dict{complexid: set([protein1, protein2,...]), .. }
    # first col: complex id
    # third col: protein id
    if format_single_protein:
        complexes = {}
        for l in ut.load_tab_file(filename):
            complexes.setdefault(l[0],set([])).add(l[2])
    else:
        complexes = dict([(l[0],set(l[1:])) for l in
                          ut.load_list_of_lists(filename)])
    return complexes
            
    
def _load_intdict(filename):
    # complexes: dict{complexid: set([protein1, protein2,...]), .. }
    # this makes a dictionary{protein1: set([protein2, protein3]), ...}
    # every interaction is found twice here for fast interaction checking
    complexes = load_complexes(filename)
    interactions = {}
    for complex,protein_set in complexes.items():
        for p in protein_set:
            partners = protein_set.copy()
            partners.remove(p)
            [interactions.setdefault(p,set([])).add(par) for par in partners]
    return interactions
