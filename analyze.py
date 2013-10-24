import compare as cp
import corum as co
import examples as ex
import features as fe
import ppi
import pairdict as pd
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

