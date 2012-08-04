from __future__ import division
import utils as ut
import pairdict as pd
from pairdict import PairDict
import ppi
import corum as co
import itertools as it

def cyto_prep(preds, cxs, arrtrain, fname, species='Hs', seqdb='ensp', 
        negmult=50):
    preds = preds_gold_standard(preds, arrtrain, species, seqdb)
    preds = preds_as_cxs(preds, cxs)
    export_ints(preds, fname, negmult)
    export_idconvert(preds, fname)

def export_ints(tested, fname, negmult, header=['id1', 'id2', 'score', 'corum', 
                'complex']):
    ut.write_tab_file([header] + [[t[0], t[1], ut.rescale(t[2],negmult)]
        + t[3:len(header)] for t in tested], fname)

def export_idconvert(preds, fname):
    pfx_convert = [[p[i], p[i].split('_')[1]] for p in preds for i in 0,1]
    pfx_convert = [['nodeid', 'ENSPID']] + pfx_convert
    ut.write_tab_file(pfx_convert, ut.pre_ext(fname,'pfx_convert'))

def preds_gold_standard(preds, arrtrain, species, seqdb):
    pdpreds = PairDict([p[:3] for p in preds])
    print len(pdpreds.d), "predicted interactions"
    ppi_cxs,_ = ppi.load_training_complexes(species, seqdb)
    pdcorum = PairDict([(i[0],i[1],'gold') for i in
                        co.pairs_from_complexes(ppi_cxs)])
    print len(pdcorum.d), "total gold standard"
    pdcomb = pd.pd_union_disjoint_vals(pdpreds, pdcorum)
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
            "preds in gold not train")
    print len([1 for p in comblist if p[2] and p[3]=='train']), "preds in train"
    # only return those that are predictions
    return [p for p in comblist if p[2]]

def preds_label_cxs(preds, cxs):
    """
    Add a new column to preds and label it with a list of arbitrary complex
    ids.
    """
    pdpreds = PairDict([p+[[]] for p in preds])
    for i,c in enumerate(cxs):
        for pair in it.combinations(c,2):
            targ = pdpreds.find(pair)
            if targ: pdpreds.d[targ][-1].append(i)
    return [list(k)+v for k,v in pdpreds.d.items()]

def preds_as_cxs(preds, cxs):
    """
    Use the complex number to both prefix the protein id and add as an
    attribute. Copy the original ids to the end for mapping in cytoscape.
    """
    preds = preds_label_cxs(preds, cxs)
    def pfx(id, cxnum):
        return str(cxnum) + '_' + id
    return [[pfx(id1,cx), pfx(id2,cx), score, corum, cx]
            for id1, id2, score, corum, cxs in preds
            for cx in cxs]
