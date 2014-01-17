import numpy as np
import compare as cp
import corum as co
import elution as el
import examples as ex
import features as fe
import orth
import ppi
import pairdict as pd
import utils as ut
from Struct import Struct

def triple_venn_consv():
    hints = co.load_havug_ints()
    ppi_cxs, clean_cxs , corconsv = ppi.load_training_complexes('Hs','Dm')
    cints = co.pairs_from_complexes(ut.i1(ppi_cxs)) # exclude huge ones
    ints23 = ut.loadpy(ut.bigd('../23_collapsenodes/Hs_filtorth025_withsc_2sp_refilt2sp_cxs_cxppis_clust27_532cxs'))[1]
    ints3 = [cp.consv_pairs(i,h2d) for i in ints23,hints,cints]
    cp.triple_venn(ints3,['map23','havug','corum'])

def correlation_enrichment(cxs):
    """
    Get the scored array for evaluating correlation of interactions within and
    between the provided set of complexes. Use this as base_exs for
    ppi.feature_array().
    """
    pdtrain, splits = ex.base_examples_single(cxs, cxs, cxs, [0,1])
    ntest_pos = len([v for k,v in pdtrain.d.items() if v[0]==1])
    base_exs = Struct(pdtrain=pdtrain, splits=splits, ntest_pos=ntest_pos)
    return base_exs

def arrfeats_prep_all_data(arrfeats, ppis, sp='Hs', gold_consv='Dm', cutoff=0.5):
    print "Adding species summary."
    arrfeats = fe.arr_add_spsummary(arrfeats, cutoff)
    print "Adding ppis."
    arrfeats = fe.arrfeats_add_ppis(arrfeats, ppis)
    _,_,all_cxs = ppi.load_training_complexes(sp, None,gold_consv) 
    pdgold = pd.PairDict(co.pairs_from_complexes(ut.i1(all_cxs)))
    print "Setting trues."
    arrfeats = fe.arrfeats_set_gold(arrfeats, pdgold) 
    return arrfeats

def enrichment_array_combined(sp_base, sp_dict_elutfs, cxs, func=np.average, 
        nsp=1, scores=['poisson']):
    """
    sp_dict_elutfs: {'Ce': [Ce_elution_1, Ce_elution_2, ...] , ...}
    """
    exs = correlation_enrichment([(i,set(c)) for i,c in enumerate(cxs)])
    elutfs = ut.flatten([elutfs for sp,elutfs in sp_dict_elutfs.items()])
    ppio = ppi.feature_array(sp_base, elutfs, exs, nsp, scores=scores, 
            extdata=[], do_filter=False)
    newarr = ppio.arrfeats
    for sp in sp_dict_elutfs.keys():
        newarr = fe.merge_features(newarr, '%s.*' % sp, func, False)
    return newarr

def pairs_found(pairs, bag):
    return len([1 for a,b in pairs if a in bag and b in bag])

def pairs_orth_found(pairs, odict, bag):
    if not odict:
        return pairs_found(pairs, bag)
    else:
        return len([1 for a,b in pairs 
            if (a in odict and b in odict 
                and len([o for o in odict[a] if o in bag]) > 0 
                and len([o for o in odict[b] if o in bag]) > 0)])

def pairs_notfound_sps(df, fs, sps='Hs Mm Sp Dm Ce'.split()):
    """
    df: dataframe with id1, id2, and the sp_evidence columns.
    fs: all the elution filenames
    """
    results = []
    for sp in sps:
        pairs = [(r['id1'],r['id2']) for i,r in
                df[df[sp+"_evidence"]!='frac'].iterrows()]
        print "%s pairs for %s" % (len(pairs), sp)
        odict = orth.odict('Hs',sp)
        orths = pairs_found(pairs, odict) if odict else len(pairs) # same sp
        fs_sp = [f for f in fs if f.find(sp+"_") > -1]
        print "%s fractionations for %s" % (len(fs_sp), sp)
        allps = el.all_prots(fs_sp)
        counts = pairs_orth_found(pairs, odict, allps)
        results.append((len(pairs), orths, counts))
    return sps, results

def barplot_foundsps(sps, counts, total_length):
    forbar = [(total_length-nofrac, nofrac-withorth, withorth-orthnoobs,orthnoobs) 
        for nofrac, withorth, orthnoobs in counts]
    import plotting as pl
    pl.stacked_bar(sps, forbar)

