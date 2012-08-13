from __future__ import division
import utils as ut
import pairdict as pd
from pairdict import PairDict
import ppi
import corum as co
import itertools as it

def cyto_prep(ppis, arrtrain, fname, f_geneatts, cxs=None, species='Hs', 
        negmult=50, pd_spcounts=None):
    header=['id1', 'id2', 'score', 'corum']
    ppis = ppis_gold_standard(ppis, arrtrain, species)
    if pd_spcounts: 
        ppis = ppis_add_counts(ppis, pd_spcounts)
        header.append('CountOfSpecies')
    if cxs:
        ppis = ppis_as_cxs(ppis, cxs)
        header.append('complexid')
    export_ints(ppis, fname, negmult, header)
    export_idconvert(ppis, cx_labels(cxs, f_geneatts), fname)


def export_ints(tested, fname, negmult, header):
    ut.write_tab_file([header] + [[t[0], t[1], ut.rescale(t[2],negmult)]
        + list(t[3:]) for t in tested], fname)

def export_idconvert(ppis, dict_cxlabels, fname):
    cxs_labeled = set([])
    pfx_convert = []
    for p in ppis:
        for i in 0,1:
            combid = p[i]
            cxid,pid = combid.split('_')
            if cxid not in cxs_labeled:
                cxlabel = dict_cxlabels[cxid]
                cxs_labeled.add(cxid)
            else:
                cxlabel = ''
            pfx_convert.append([combid, pid, cxid, cxlabel])
    pfx_convert = [['nodeid', 'ENSGID', 'complexid', 'ComplexLabel']] \
            + pfx_convert
    ut.write_tab_file(pfx_convert, ut.pre_ext(fname,'pfx_convert'))

def cx_labels(cxs, f_geneatts):
    """
    cxs is just a list of lists of gene ids.
    f_geneatts points to a text file with col0=geneid, col1=name, 
    col2=description.
    Return a dict of {cxid: label}
    """
    id2name = dict([(p[0],p[1]) for p in ut.load_tab_file(f_geneatts)])
    d = {}
    print len(set([g for genes in cxs for g in genes if g not in id2name])), \
            "genes without name"
    for cxid,genes in enumerate(cxs):
        gnames = [id2name[g] for g in genes if g in id2name]
        matches = matches_from_names(gnames)
        usen = 3
        topn = matches[:usen]
        ntotal = len(genes)
        nother = ntotal - sum([m[0] for m in topn])
        d[str(cxid)] = ' '.join([match for n,match in topn]) + (' +%s' % nother 
                if nother > 0 else '')
    return d

def matches_from_names(innames):
    names = sorted(innames)
    i = 0
    matches = []
    while i < len(names):
        nmatches, match = best_match(names,i)
        i += nmatches
        matches.append((nmatches,match))
    return sorted(matches, key=lambda x:x[0], reverse=True)

def best_match(names, i):
    """
    Find matches starting from index i.
    Matches have to be at least 3 letters, otherwise it's a mess for big sets.
    Go with largest set of matches instead of longer matches.
    Return (number of matches, matching chars)
    """
    matchlen = 3
    nmatches = 1
    while (i+nmatches < len(names) and names[i+nmatches][:matchlen] ==
            names[i][:matchlen]):
        nmatches += 1
    while matchlen < 8 and all([names[i+j][:matchlen+1] ==
        names[i][:matchlen+1] for j in range(1,nmatches)]):
        matchlen += 1
    return nmatches, names[i][:matchlen]

def ppis_gold_standard(ppis, arrtrain, species):
    pdppis = PairDict([p[:3] for p in ppis])
    print len(pdppis.d), "predicted interactions"
    ppi_cxs,_ = ppi.load_training_complexes(species)
    pdcorum = PairDict([(i[0],i[1],'gold') for i in
                        co.pairs_from_complexes(ppi_cxs)])
    print len(pdcorum.d), "total gold standard"
    pdcomb = pd.pd_union_disjoint_vals(pdppis, pdcorum)
    pdtrainpos = PairDict([(t[0],t[1]) for t in arrtrain if t[2]==1])
    print len(pdtrainpos.d), "total train interactions"
    counterrs = 0
    for tpair in pdtrainpos.d:
        cpair = pdcomb.find(tpair)
        if pdcomb.d[cpair][1] != 'gold':
            #print 'error: train should be subset', tpair
            counterrs += 1
        else:
            pdcomb.d[cpair][1] = 'train'
    if counterrs: print "number of training not found in gold std:", counterrs
    comblist = [list(k)+list(v) for k,v in pdcomb.d.items()]
    print (len([1 for p in comblist if p[2] and p[3]=='gold']), 
            "ppis in gold not train")
    print len([1 for p in comblist if p[2] and p[3]=='train']), "ppis in train"
    # only return those that are predictions
    return [p for p in comblist if p[2]]

def ppis_add_counts(ppis, pd_spcounts):
    n_notfound = 0
    newppis = []
    for p in ppis:
        pair = pd_spcounts.find((p[0], p[1]))
        if not pair: n_notfound += 1
        newppis.append(p + pd_spcounts.d.get(pair,[-1])) # already a list
    print "appending species counts;", n_notfound, "not found."
    return newppis

def ppis_label_cxs(ppis, cxs):
    """
    Add a new column to ppis and label it with a list of made-up complex
    ids.
    """
    pdppis = PairDict([p+[[]] for p in ppis])
    for i,c in enumerate(cxs):
        for pair in it.combinations(c,2):
            targ = pdppis.find(pair)
            if targ: pdppis.d[targ][-1].append(i)
    return [list(k)+v for k,v in pdppis.d.items()]

def ppis_as_cxs(ppis, cxs):
    """
    Use the complex number to both prefix the protein id and add as an
    attribute. Copy the original ids to the end for mapping in cytoscape.
    """
    ppis = ppis_label_cxs(ppis, cxs)
    # Requires that the last column be a list of complex ids. Replace that.
    def pfx(id, cxnum):
        return str(cxnum) + '_' + id
    return [[pfx(p[0],cx), pfx(p[1],cx)] + p[2:-1] + [cx]
            for p in ppis for cx in p[-1]]

