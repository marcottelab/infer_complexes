from __future__ import division
import numpy as np
import random
import os
import utils as ut
import ppi
import seqs

def pairs(fname):
    return [list(e) for e in ut.load_tab_file(fname)]

def load_corum(fname):
    """
    Returns a list of tuples: (name, set(uniprotIDs), species)
    """
    return [(l[1], set(l[4].split(',')), l[3]) 
            for l in ut.load_list_of_lists(fname, sep=';')[1:]]

def load_havug_ppis():
    hints = ut.load_list_of_lists('../../docs/SupplementaryTableS2.tab')
    u2e = ut.dict_inverse_sets(ut.load_dict_sets('../../data/convert/Hs2Hs_uni.tab'))
    hints = [[list(u2e.get(p,['NoTranslation']))[0] for p in c[:2]]+[c[2]] for c in hints]
    hcxsu = ut.load_list_of_type('havig_complexes.tab',set)
    hcxs = ut.i1(co.convert_complexes(dict([(i,c) for i,c in enumerate(hcxsu)]), u2e, seqs.load_prots_from_fasta('../../data/sequences/canon/Hs.fasta')))
    hints = cl._filter_ints(hints, hcxs)
    return hints


def load_ppi_cxs(minlen=2, maxlen=50, sp_match='Human'):
    """
    Returns a list of sets of uniprot ids.
    No de-duplication or anything else.
    Expected that this list may have duplicates.
    """
    fname = ut.proj_path('corum_cxs')
    return [(name,ps) for name,ps,s in load_corum(fname) 
            if (not sp_match or s==sp_match) 
            and len(ps)>=minlen and len(ps)<=maxlen]

def load_merged_cxs(cutoff=0.5, func=None, sep='|', **kwargs):
    if func==None:
        import compare as cp
        func = cp.simpson
    cxs = load_ppi_cxs(**kwargs)
    return merge_atonce(cxs, cutoff, func, sep)
    #if atonce: # deprecated since I didn't fix merge_iter for handling names
    #else:
        #return merge_iter(cxs, cutoff, func)

def merge_atonce(psets, cutoff, func, sep):
    """
    Once difference seems to be that say a series of overlapping doubles can't
    string together with this approach, whereas it can with the iter approach.
    """
    to_merge = list(psets)
    merged = []
    while len(to_merge)>1:
        c1,c1ps = random.choice(to_merge)
        to_merge.remove((c1,c1ps))
        matches = [(c,ps) for c,ps in to_merge if func(ps,c1ps)>cutoff]
        for m in matches: to_merge.remove(m)
        newname = sep.join([c1]+ut.i0(matches))
        newps = reduce(set.union, [c1ps]+ut.i1(matches))
        merged.append((newname,newps))
    merged.append(to_merge[0])
    return merged

# not fixed for handling names
#def merge_iter(psets, cutoff, func):
    #merged = list(psets)
    #repeat = True
    #while repeat: # iterate until no more merging can happen
        #repeat = False
        #to_merge = merged
        #merged = []
        #while len(to_merge)>1: # single loop through
            #c1 = random.choice(to_merge)
            #to_merge.remove(c1)
            #c2 = cp.best_match_item(c1, to_merge, func)
            #if func(c1,c2) > cutoff:
                #merged.append(set.union(c1,c2))
                #to_merge.remove(c2)
                #repeat = True
            #else:
                #merged.append(c1)
        #merged.append(to_merge[0]) # base case -- last item
    #return merged


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
    Convert a list of complexes from one proteinid scheme to another.
    convert: dict of { fromid1: set([toid1, toid2, ...]), ...}
    only_to_prots: a _set_ of outids, such that we only convert to any in that set
    If only_to_prots is None, we use all the outids.
    For brevity in code I assumed wlog from uniprot 'u' to ensp 'e'.
    """
    complexes = dict(complexes)
    assert type(only_to_prots) == type(set([])), 'Prots must be set for speed'
    convert_filtered = dict([(u,[e for e in convert[u] if (only_to_prots is
        None or e in only_to_prots)][0]) for u in convert if len([e for e in
        convert[u] if (only_to_prots is None or e in only_to_prots)])>0])
    out_complexes = dict([(c,set([convert_filtered[u] for u in complexes[c] if
        u in convert_filtered])) for c in complexes if len([convert_filtered[u]
        for u in complexes[c] if u in convert_filtered])>1])
    return out_complexes.items()

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

class CLookup(object):

    def __init__(self, cxs=None, consv_sp=''):
        gtrans = seqs.GTrans()
        if cxs is None:
            cxs,_,_ = ppi.load_training_complexes('Hs',consv_sp)
        allps = reduce(set.union, [set(c[1]) for c in cxs])
        notfound = len([p for p in allps if p not in gtrans.id2name])
        if notfound > 0: print '%s gene names not found' % notfound
        self.cxs = [(c[0], set([gtrans.id2name.get(p,'noname') for p in
            c[1]])) for c in cxs]

    def findcs(self,gname):
        return [p for p in self.cxs if gname.lower() in p[1]]
