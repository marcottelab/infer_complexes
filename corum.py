from __future__ import division
import numpy as np
import random
import os
import utils as ut

keys = {
    #'uppi': ('~/Dropbox/complex/data/corum/corum_overlaps_forppi.tab', 'lists'),
    #'eppi': ('~/Dropbox/complex/data/corum/corum_overlaps_forppi_en.tab', 'lists'),
    #'uclean': ('~/Dropbox/complex/data/corum/corum_clean_nooverlaps_supptable3.tab', 'singles'),
    #'eclean':
    #('~/Dropbox/complex/data/corum/corum_clean_nooverlaps_supptable3_ensp.tab',
    #'singles')
    'hs_eptrain': '~/Dropbox/complex/data/corum/pairs/ensp_pairs_ppi.tab',
    'hs_eptrain_negs': '~/Dropbox/complex/data/corum/pairs/ensp_pairs_ppi_negs.tab',
    'hs_eptest':
    '~/Dropbox/complex/data/corum/pairs/ensp_pairs_ppi_exclude.tab',
    'hs_eptest': '~/Dropbox/complex/data/corum/pairs/ensp_pairs_ppi_negs_exclude.tab',
    'hs_uptrain': '~/Dropbox/complex/data/corum/pairs/uni_pairs_ppi.tab',
    'hs_uptrain_negs': '~/Dropbox/complex/data/corum/pairs/uni_pairs_ppi_negs.tab',
    'hs_uptest':
    '~/Dropbox/complex/data/corum/pairs/uni_pairs_ppi_exclude.tab',
    'hs_uptest_negs': '~/Dropbox/complex/data/corum/pairs/uni_pairs_ppi_negs_exclude.tab',
    'ce_eptrain': '~/Dropbox/complex/data/corum/pairs/ce_ensp_pairs_ppi.tab',
    'ce_eptrain_negs': '~/Dropbox/complex/data/corum/pairs/ce_ensp_pairs_ppi_negs.tab',
    'ce_eptest':
    '~/Dropbox/complex/data/corum/pairs/ce_ensp_pairs_ppi_exclude.tab',
    'ce_eptest_negs': '~/Dropbox/complex/data/corum/pairs/ce_ensp_pairs_ppi_negs_exclude.tab'
    }


def _dep_pairs_key(key):
    return [list(e) for e in ut.load_tab_file(os.path.expanduser(keys[key]))]

def pairs(fname):
    return [list(e) for e in ut.load_tab_file(fname)]


def pairs_from_complexes(complexes):
    intdict = corum_ints_duped(complexes)
    int_dedup = ut.dict_dedup(intdict)
    pairs = []
    for p, partners in int_dedup.items():
        for par in partners:
            pairs.append((p,par))
    return pairs

def load_complexes_singleline(filename, startindex=1):
    """
    (Usually for 'ppi' overlapping complexes)
    """
    # load corum-type file into a dictionary
    # complexes: dict{complexid: set([protein1, protein2,...]), .. }
    # first col: complex id
    filename = os.path.expanduser(filename)
    complexes = {}
    # PPI complex set from traver has duplicate complex names with different
    # members.  Using this approach means all the members from any lines with
    # that complex's name get added.
    for l in ut.load_tab_file(filename):
        for i in l[startindex:]:
            complexes.setdefault(l[0],set([])).add(i)
    return complexes

def load_complexes_multiline(filename):
    """
    (Usually for 'clean' complexes).
    Load complexes in a file in the style of supp table 3: complexid,
    complexname, singlemember.
    """
    filename = os.path.expanduser(filename)
    complexes = {}
    for l in ut.load_tab_file(filename):
        complexes.setdefault(l[1],set([])).add(l[2])
    return complexes
            
    
def corum_ints_duped(complexes):
    # complexes: dict{complexid: set([protein1, protein2,...]), .. }
    # this makes a dictionary{protein1: set([protein2, protein3]), ...}
    # every interaction is found twice here for fast interaction checking
    interactions = {}
    for complex,protein_set in complexes.items():
        for p in protein_set:
            partners = protein_set.copy()
            partners.remove(p)
            [interactions.setdefault(p,set([])).add(par) for par in partners]
    return interactions

def divide_pairs(pairs, fractions):
    # fractions like [.35,.5,.85,1]
    cv_sets = [[] for f in fractions]
    for p in pairs:
        r = random.random()
        cv_sets[ np.min(np.where(r < np.array(fractions))) ].append(p)
    return cv_sets


def convert_complexes(complexes, convert, only_to_prots=None):
    """
    Convert a dict of complexes from one proteinid scheme to another.
    convert: dict of { fromid1: set([toid1, toid2, ...]), ...}
    only_to_prots: a _set_ of outids, such that we only convert to any in that set
    If only_to_prots is None, we use all the outids.
    For brevity in code I assumed wlog from uniprot 'u' to ensp 'e'.
    """
    assert type(only_to_prots) == type(set([])), 'Prots must be set for speed'
    convert_filtered = dict([(u,[e for e in convert[u] if (only_to_prots is
        None or e in only_to_prots)][0]) for u in convert if len([e for e in
        convert[u] if (only_to_prots is None or e in only_to_prots)])>0])
    out_complexes = dict([(c,set([convert_filtered[u] for u in complexes[c] if
        u in convert_filtered])) for c in complexes if len([convert_filtered[u]
        for u in complexes[c] if u in convert_filtered])>1])
    return out_complexes

def write_pos_neg_pairs(complexes, complexes_exclude, fname_pos):
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
