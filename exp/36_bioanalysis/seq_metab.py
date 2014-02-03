from __future__ import division
import random
from scipy import io
from scipy.stats import ks_2samp
import numpy as np
import re
import orth
import pairdict as pd
import ppi_utils as pu
import utils as ut

def figures(recon_fname, exclude_recon_fname, kegg_fname, all_ppis_fname, 
        recon_pairs=None):
    rpairs = recon_pairs or load_seq_pairs(recon_fname,
            ut.load_list(exclude_recon_fname))
    kpairs = load_kegg_sequentials(kegg_fname)
    pdk = pd.PairDict(kpairs)
    intpairs = [p for p in rpairs if pdk.contains(p)]
    ppis = pu.load_ppis(all_ppis_fname) 
    plot_pairs_randoms_etc(intpairs, ppis)
    plot_cdf_pos_randoms(intpairs, ppis)

def plot_cdf_pos_randoms(pospairs, ppis):
    import plotting as pl
    pl.figure()
    pos,neg1 = pl.hist_pairs_nonpairs(ppis, pospairs, negmult=1, do_plot=False)
    pos,neg100 = pl.hist_pairs_nonpairs(ppis, pospairs, negmult=100,
            do_plot=False)
    for pairs in pos, neg1, neg100:
        pl.cdf(pairs,bins=np.arange(0,1.01,.01))
    pl.xlabel("PPI score")
    pl.ylabel("Cumulative fraction of population")
    pl.title('Several percent of sequential enzymes are high-scoring,\ncompared to much less than one percent for random shufflings')
    pl.legend(['Sequentials','Size-matched reshuffled','100x larger set of reshuffled'],loc=4)


def load_kegg_brite(fname):
    """
    The pathways dict is keyed by pathway names, and each value is a list of
    enzyme pairs, which is a set of enzyme ids and a set of gene ids.
    """
    lines = [L.strip().split('\t')[0] 
            for L in file(fname, 'r') if re.match("^[A-Z]", L)]
    d_pathways = {}
    d_genes = {}
    pathway = 'Not Found'
    for L in lines:
        #if L[0]=='B':
            #category = L[L.find('<b>')+3:L.find('</b>')]
        contents = L[1:].strip()
        if L[0]=='C':
            pathway = contents
            d_pathways[pathway] = []
            prior_enzymes = set([])
        elif L[0]=='D':
            entrez = contents.split()[0]
            if contents.find('(EC:') == -1:
                enzymes = set(['Unknown'])
            else:
                enzymes = set(contents[contents.find('(EC:')+4:-1].split())
            if len(set.intersection(enzymes, prior_enzymes)) == 0:
                d_pathways[pathway].append([enzymes, set([])])
                prior_enzymes = enzymes
            else:
                d_pathways[pathway][-1][0] = set.union(enzymes, prior_enzymes)
                prior_enzymes = set.union(enzymes, prior_enzymes)
            d_pathways[pathway][-1][1].add(entrez)
    return d_pathways

def load_kegg_sequentials(fname, do_convert=True):
    dkegg = load_kegg_brite(fname)
    kegg_paths = [ut.i1(v) for v in dkegg.values() if v]
    def path_pairs(list_path):
        return [(list_path[i],list_path[i+1]) for i in range(len(list_path)-1)]
    group_pairs = ut.flatten([path_pairs(lpath) for lpath in kegg_paths])
    #if return_groups:
        #if conv_dict:
            #return convert_groups_singles(labeled_pairs, conv_dict)
        #else:
            #return labeled_pairs
    single_pairs = [(xi,yi) for x,y in group_pairs for xi in x for yi in y]
    unique_pairs = pu.dedupe(single_pairs)
    print "%s total, %s single, %s unique pairs returned" % (
            len(group_pairs), len(single_pairs), len(unique_pairs)) 
    if do_convert:
        conv_dict = ut.dict_inverse_sets(orth.convert_dict('Hs','Hs_entrez'))
        conv_pairs = convert_pairs_singles(unique_pairs, conv_dict)
        print "%s converted pairs with 1-1 matches" % len(conv_pairs)
        return conv_pairs
    else:
        return unique_pairs


def random_pairs(pairs, npairs):
    ps = list(set(ut.i0(pairs) + ut.i1(pairs)))
    rpairs = [(random.choice(ps), random.choice(ps)) for i in
            range(int(npairs*1.5))]
    return pu.dedupe(rpairs)[:npairs]

def entrez_desc():
    return dict([(l[4][2:],l[2].split('[')[0]) for l in
        ut.load_lol(ut.config('gene_desc_Hs'))])

def adj_matrix(S):
    St = S.transpose()
    St = St > 0
    S = S < 0
    mat_adj = np.matrix(St)*np.matrix(S)
    return mat_adj
#def adj_matrix(S):
    #St = S.transpose()
    #St = St.todense() > 0
    #S = S.todense() < 0
    #mat_adj = St*S
    #return mat_adj

def seq_pairs(S, enzymes, mat_adj=None, conv_dict=None, return_groups=False):
    """
    - S: stoichiometric matrix of reactants (-) and products (+)
    - enzymes: list of enzymes corresponding to columns in S
    """
    mat_adj = mat_adj if mat_adj is not None else adj_matrix(S)
    rows, cols = np.where(mat_adj > 0)
    #wrestling with where output
    rows, cols = [np.array(x)[0] for x in rows, cols] 
    rowcols = ut.zip_exact(rows, cols)
    # dedup, keeping only upper off-diagonal
    rowcols = [(row,col) for row,col in rowcols if row < col]
    pairs = [(enzymes[row],enzymes[col]) for row,col in rowcols]
    # filter out blanks
    labeled_pairs = [(x,y) for x,y in pairs if x and y]
    if return_groups:
        if conv_dict:
            return convert_groups_singles(labeled_pairs, conv_dict)
        else:
            return labeled_pairs
    single_pairs = [(xi,yi) for x,y in labeled_pairs 
            for xi in x.split() for yi in y.split()]
    unique_pairs = pu.dedupe(single_pairs)
    print "%s total, %s labeled, %s single, %s unique pairs returned" % (len(pairs), 
            len(labeled_pairs), len(single_pairs), len(unique_pairs)) 
    if conv_dict:
        conv_pairs = convert_pairs_singles(unique_pairs, conv_dict)
        print "%s converted pairs with 1-1 matches" % len(conv_pairs)
        return conv_pairs
    else:
        return unique_pairs

def grRules_to_enzymes(grRules):
    return [' '.join([x.strip('()').split('.')[0] for x in c.split(' or ')]) 
            for c in grRules]

def convert_groups_singles(pairgroups, cdict):
    def conv_group(items, cdict):
        return [x for x in (conv_single(i, cdict) for i in items) if x]
    return [(a,b) for a,b in 
            ((conv_group(x, cdict), conv_group(y, cdict)) for x,y in pairgroups) 
            if len(a) > 0 and len(b) > 0]

def convert_pairs_singles(pairs, cdict):
    """
    Return pairs where both items are found and have only single matches in the
    provided conversion dict, mapping items to sets.
    """
    return [(a,b) for a,b in 
            ((conv_single(x, cdict),conv_single(y, cdict)) for x,y in pairs) 
            if a and b]

def conv_single(a, conv_dict_sets):
    b = conv_dict_sets.get(a,set([]))
    return list(b)[0] if len(b) == 1 else None

def load_metabolic_data(fname):
    def unnest(nestarr):
        return [x[0][0] if x else '' for x in nestarr]
    mobjs = io.loadmat(fname)
    S = np.asarray(mobjs['S'].todense())
    entrez_enzymes = grRules_to_enzymes(unnest(mobjs['grRules']))
    rxn_names = unnest(mobjs['rxnNames'])
    met_names = unnest(mobjs['metNames'])
    return S, entrez_enzymes, rxn_names, met_names

def load_seq_pairs(fname, metab_exclude=None):
    """
    metab_exclude: should be in sequential_metab/metabolites_exclude.txt
    """
    S, entrez_enzymes, rnames, mnames = load_metabolic_data(fname)
    ez2en = ut.dict_inverse_sets(orth.convert_dict('Hs','Hs_entrez'))
    if metab_exclude:
        print "Excluding %s metabolites, filtering rxns" % len(metab_exclude)
        S, entrez_enzymes = filter_rxns_metabs(S, entrez_enzymes, rnames,
                mnames, metab_exclude) 
    else:
        print "No filtering of metabolites and rxns."
    sequentials = seq_pairs(S, entrez_enzymes, conv_dict=ez2en)
    return sequentials

def filter_rxns_metabs(S, enzymes, rnames, mnames, metab_exclude):
    rinds = filtered_rxn_inds(rnames)
    metab_exclude = set(metab_exclude)
    metab_exclude.add('')
    minds = [i for i,m in enumerate(mnames) if m not in metab_exclude]
    print "Keeping %s of %s reactions" % (len(rinds), len(rnames))
    print "Keeping %s of %s metabolites" % (len(minds), len(mnames))
    newS = S[minds,:]
    newS = newS[:,rinds]
    newenz = ut.list_inds(enzymes, rinds)
    return newS, newenz
    

def plusn(sequentials, limit='auto'):
    """
    From pairs of sequentials, return pairs of items separated by n steps.
    """
    def hop_across(pairsa, pairsb, pairs_exclude):
        a,b = [ut.dict_sets_from_tuples(p) for p in pairsa,pairsb]
        pd_exclude = pd.PairDict(pairs_exclude)
        newpairs = [(x,z) for x in a for y in a.get(x,[]) for z in b.get(y,[]) 
                if not (x==z or pd_exclude.contains((x,z)))]
        newpairs = pu.dedupe(newpairs)
        return newpairs
    if limit=='auto':
        limit = len(sequentials)
    def maybe_sample(items, limit):
        return random.sample(items, limit) if len(items)>limit else items
    plus2s = hop_across(sequentials, sequentials, sequentials)
    #limited_plus2s = random.sample(plus2s, limit) if len(plus2s)>limit else plus2s
    print len(plus2s), "n+2s found"
    plus3s = hop_across(sequentials, plus2s, sequentials + plus2s)
    #plus3s = hop_across(sequentials, limited_plus2s, sequentials + plus2s)
    print len(plus3s), "n+3s found"
    #limited_plus3s = random.sample(plus3s, limit) if len(plus3s)>limit else plus3s
    #plus3s = [(x,z) for x in adj_dict for y in adj_dict.get(x,[]) for z in
            #dict_plus2s.get(y,[]) if not pd_exclude.contains((x,z))]
    plus4s = hop_across(plus2s, plus2s, sequentials + plus2s +
            plus3s)
    #plus4s = hop_across(limited_plus2s, limited_plus2s, sequentials + plus2s +
            #plus3s)
    print len(plus4s), "n+4s found"
    #limited_plus4s = random.sample(plus4s, limit) if len(plus4s)>limit else plus4s
    #plus4s = [(x,z) for x in dict_plus2s for y in dict_plus2s.get(x,[]) for z in
            #dict_plus2s.get(y,[]) if not pd_exclude.contains((x,z))]
    return [maybe_sample(items, limit) for items in plus2s, plus3s+plus4s]

def plus2(sequentials, limit='auto'):
    """
    From pairs of sequentials, return pairs of items separated by n steps.
    """
    def hop_across(pairsa, pairsb, pairs_exclude):
        a,b = [ut.dict_sets_from_tuples(p) for p in pairsa,pairsb]
        pd_exclude = pd.PairDict(pairs_exclude)
        newpairs = [(x,z) for x in a for y in a.get(x,[]) for z in b.get(y,[]) 
                if not (x==z or pd_exclude.contains((x,z)))]
        newpairs = pu.dedupe(newpairs)
        return newpairs
    if limit=='auto':
        limit = len(sequentials)
    def maybe_sample(items, limit):
        return random.sample(items, limit) if len(items)>limit else items
    plus2s = hop_across(sequentials, sequentials, sequentials)
    #limited_plus2s = random.sample(plus2s, limit) if len(plus2s)>limit else plus2s
    print len(plus2s), "n+2s found"
    #plus3s = hop_across(sequentials, plus2s, sequentials + plus2s)
    #plus3s = hop_across(sequentials, limited_plus2s, sequentials + plus2s)
    #print len(plus3s), "n+3s found"
    #limited_plus3s = random.sample(plus3s, limit) if len(plus3s)>limit else plus3s
    #plus3s = [(x,z) for x in adj_dict for y in adj_dict.get(x,[]) for z in
            #dict_plus2s.get(y,[]) if not pd_exclude.contains((x,z))]
    #plus4s = hop_across(plus2s, plus2s, sequentials + plus2s +
            #plus3s)
    #plus4s = hop_across(limited_plus2s, limited_plus2s, sequentials + plus2s +
            #plus3s)
    #print len(plus4s), "n+4s found"
    #limited_plus4s = random.sample(plus4s, limit) if len(plus4s)>limit else plus4s
    #plus4s = [(x,z) for x in dict_plus2s for y in dict_plus2s.get(x,[]) for z in
            #dict_plus2s.get(y,[]) if not pd_exclude.contains((x,z))]
    #return [maybe_sample(items, limit) for items in plus2s, plus3s, plus4s]
    return maybe_sample(plus2s, limit)

def plot_pairs_randoms_etc(sequentials, score_ppis, plusns=None):
    import plotting as pl
    pl.figure()
    plus2s, plus3s4s = plusns or plusn(sequentials)
    #plus2s = plusns or plusn(sequentials)
    rand_pairs = random_pairs(sequentials, len(sequentials))
    scores = [pl.hist_pairs_nonpairs(score_ppis, pairs, negmult=10,
        do_plot=False)
        for pairs in sequentials, plus2s, plus3s4s, rand_pairs]
        #for pairs in sequentials, plus2s, plus3s, plus4s, rand_pairs]
    ks_pvals = [ks_2samp(pos,neg)[1] for pos,neg in scores] # [1] is p-value
    logvals = [-np.log10(pval) for pval in ks_pvals]
    #pl.bar_plot(['%s\np < %0.3g' % (x,y) for x,y in zip('n,n+1 n,n+2 n,n+3 random'.split(), ks_pvals)], logvals)
    pl.bar_plot(['%s\np < %0.3g' % (x,y) for x,y in zip('n,n+1 n,n+2 n,n+3/4 random'.split(), ks_pvals)], logvals)
    pl.ylabel('-log10(p-value) : two-sample K/S test')
    pl.title('Intersection of recon and kegg sequential pairs\nNpairs=%s; %s n+2s, %s n+3s,n+4s' % (len(sequentials), len(plus2s), len(plus3s4s)))
    return logvals

def common_metabolites(fname, rxn_inds=None, printn = 30):
     mobjs = io.loadmat(fname) 
     S = mobjs['S']
     Sdense = S.todense()
     Sarr = np.asarray(Sdense)
     if rxn_inds:
         Sarr = Sarr[:, rxn_inds]
     Sbin = Sarr != 0
     binsums = Sbin.sum(axis=1)
     msums = [(m[0][0] if m else '-',mf[0][0] if mf else '-',s) for m,mf,s in zip(mobjs['metNames'],mobjs['metFormulas'],binsums)]
     msums.sort(key=lambda x: x[2], reverse=True)
     ut.print_lol(msums[:printn])
     return msums

#def exclude_metabolites(fname=None, metab_sums=None, limit):
     #metab_sums = metab_sums or common_metabolites(fname)

def filter_S_rxns(S, exclude_list):
    pass

exclude_rxns_default = ['transport', 'diffusion', 'degradation', 'flux', 
        'protein assembly']

def filtered_rxn_inds(rxnNames, exclude_list=exclude_rxns_default):
    """
    Exclude if the name.lower() contains anything in the exclude list.
    Return the indices to keep.
    """
    return [i for i,name in enumerate(rxnNames) 
            if max([name.lower().find(ex) for ex in exclude_list]) == -1]

