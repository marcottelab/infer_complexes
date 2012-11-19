import compare as cp
import corum as co

def triple_venn_consv():
    hints = co.load_havug_ints()
    ppi_cxs, clean_cxs , corconsv = ppi.load_training_complexes('Hs','Dm')
    # cints = co.pairs_from_complexes(dict(corconsv))
    cints = co.pairs_from_complexes(dict(ppi_cxs)) # exclude huge ones
    ints23 = ut.loadpy(ut.bigd('../23_collapsenodes/Hs_filtorth025_withsc_2sp_refilt2sp_cxs_cxppis_clust27_532cxs'))[1]
    ints3 = [cp.consv_pairs(i,h2d) for i in ints23,hints,cints]
    cp.triple_venn(ints3,['map23','havug','corum'])
