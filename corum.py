from __future__ import division
from collections import defaultdict
import numpy as np
import random
import os
import utils as ut
import ppi
import seqs
import cluster as cl
from Struct import Struct


def pairs(fname):
    return [list(e) for e in ut.load_tab_file(fname)]

def load_corum(fname, filter_methods, do_dedupe):
    """
    Returns a list of tuples: (name, set(uniprotIDs), species)
    """
    lines = [l[:7] for l in ut.load_tab_file(fname, sep=';')][1:]
    cxs = [(name, set(prots.split(',')), species, method) 
            for _,name,_,species,prots,_,method in lines]
    if filter_methods:
        print "Filtering corum methods."
        keep_methods = set([x[0] for x in
            (ut.load_tab_file(ut.proj_path('corum_methods'))) if int(x[3])==1])
        cxs = [(n,p,s) for n,p,s,methods in cxs 
                if (len([m for m in methods.split('|') 
                    if m.split('-')[0].strip() in keep_methods]) > 0)]
    else:
        cxs = [(n,p,s) for n,p,s,m in cxs]
    return cxs

def load_havug_ppis():
    hints = ut.load_list_of_lists('../../docs/SupplementaryTableS2.tab')
    u2e = ut.dict_inverse_sets(ut.load_dict_sets('../../data/convert/Hs2Hs_uni.tab'))
    hints = [[list(u2e.get(p,['NoTranslation']))[0] for p in c[:2]]+[c[2]] for c in hints]
    return hints

def load_havug_cxppis():
    hints = load_havug_ppis()
    hcxs = load_havug_cxs()
    hints = cl._filter_ints(hints, ut.i1(hcxs))
    return hints

def load_havug_cxs(convert_ensg=True):
    fname = ut.proj_path('havug_cxs')
    u2e = ut.dict_inverse_sets(ut.load_dict_sets(
        '../../data/convert/Hs2Hs_uni.tab'))
    hcxs = ut.load_list_of_type(fname,set)
    if convert_ensg:
        hcxs = convert_complexes(dict([(i,c) for i,c in
            enumerate(hcxs)]), u2e,
            seqs.load_prots_from_fasta('../../data/sequences/canon/Hs.fasta'))
    return hcxs


def filter_location(cxs, go_location):
    """
    Return only those cxs for x (cyto/nuc) where the go cell compartment is either:
    - more proteins annotated with x than y
    - more than half of the proteins annotated with x
    """
    assert go_location in ['cyto','nuc'], 'location not supported'
    keys = ['cyto_prots','nuc_prots']
    yes_key,no_key = keys if go_location=='cyto' else keys[::-1]
    yes_prots,no_prots = [go_assoc_prots(ut.proj_path(key)) 
            for key in yes_key, no_key]
    return [c for c in cxs 
            if len(set.intersection(c[1],yes_prots))/len(c[1]) > .5 
            or (len(set.intersection(c[1],yes_prots)) -
                len(set.intersection(c[1],no_prots))) > 0]

def print_rib_count(cxs, label):
    print ('ribosome cxs ps %s' % label,  len([x for x in cxs if
        x[0].find('ibosom')>-1]), len(cxs), sum([len(c[1]) for c in cxs]))


def load_ppi_cxs(minlen=2, maxlen=50, sp_match='Human', go_location=None,
        do_filter_methods=True, dedupe_names=True):
    """
    Returns a list of sets of uniprot ids.
    No de-duplication or anything else.
    Expected that this list may have duplicates.
    """
    fname = ut.proj_path('corum_cxs')
    cxs = load_corum(fname, do_filter_methods, dedupe_names)
    print_rib_count(cxs, 'a')
    if sp_match: cxs = ut.list_filter_value(cxs, 2, sp_match)
    #print_rib_count(cxs, 'b')
    cxs = [c for c in cxs if len(c[1])>=minlen and len(c[1])<=maxlen]
    #print_rib_count(cxs, 'c')
    cxs = [(name,ps) for name,ps,spec in cxs]
    #print len(cxs)
    if go_location:
        print "Filtering corum by go location"
        cxs = filter_location(cxs, go_location)
    return cxs

def load_merged_cxs(cutoff=0.5, func=None, sep='|', go_location=None, **kwargs):
    if func==None:
        import compare as cp
        func = cp.simpson
    cxs = load_ppi_cxs(**kwargs)
    cxs = merge_atonce(cxs, cutoff, func, sep)
    if go_location:
        cxs = filter_location(cxs, go_location)
    return cxs
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


def convert_complexes(complexes, convert, include_set=None):
    """
    Convert a list of complexes from one proteinid scheme to another.
    convert: dict of { fromid1: set([toid1, toid2, ...]), ...}
    include_set: a _set_ of outids, such that we only convert to any in that set
    If only_to_prots is None, we use all the outids.
    For brevity in code I assumed wlog from uniprot 'u' to ensp 'e'.
    """
    #complexes = dict(complexes)
    assert type(include_set) == type(set([])), 'Prots must be set for speed'
    #convert_filtered = [(u,[e for e in convert[u] if (only_to_prots is
        #None or e in only_to_prots)][0]) for u in convert] if len([e for e in
        #convert[u] if (only_to_prots is None or e in only_to_prots)])>0])
    convert_filtered = {}
    for u,es in convert.items():
        es_filt = [e for e in es if e in include_set or include_set is None]
        if len(es_filt)>0: 
            convert_filtered[u] = es_filt[0]
    out_complexes = []
    for c,us in complexes:
        es = set([convert_filtered[u] for u in us if u in convert_filtered])
        if len(es)>0:
            out_complexes.append((c,es))
    #out_complexes = dict([(c,set([convert_filtered[u] for u in complexes[c] if
        #u in convert_filtered])) for c in complexes if len([convert_filtered[u]
        #for u in complexes[c] if u in convert_filtered])>1])
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

def go_assoc_prots(fname):
    """
    Just return the second item in each line if the line isn't commented (!).
    Assume this file is already filtered only to contain the desired set of
    proteins.
    """
    return set([l.split('\t')[1] for l in file(fname,'r') if l[0]!='!'])

def count_ints(cxs):
    """
    Given cxs, a list of [(name, set([id1,id2...])), ...], return the number of
    possible positive pairwise interactions implied.
    """
    def count_pairs(items):
        n = len(items)
        return int(n*(n-1)/2)
    return sum([count_pairs(c[1]) for c in cxs])

def keep_longest(cxs):
    d_cxs = {}
    for name, ps in cxs:
        if (name not in d_cxs) or (len(ps) > len(d_cxs[name])):
            d_cxs[name] = ps
    return d_cxs.items()

#def dedupe_names(cxs):
    #"""
    #Remove possibility of name bashing by renaming complexes that would bash.
    #NOT WORKING 20130806.
    #"""
    #names_counts = defaultdict(int)
    #print len(cxs)
    #cxs_deduped = []
    #d_cxs = {}
    #for name, ps in cxs:
        #if names_counts[name] > 0:
            #if d_cxs[name] != ps:
                #use_name = '%s_%s' % (name, names_counts[name]+1)
                #names_counts[name] += 1
            #else:
                ## skip: duplicate
                #continue
        #else:
            #use_name = name
        #d_cxs[name] = ps
        #cxs_deduped.append((use_name, ps))
    #print len(cxs_deduped)
    #return cxs_deduped

