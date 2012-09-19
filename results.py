import numpy as np
import utils as ut
from Struct import Struct
import ml
import features as fe
import cluster as cl
import ppi
import cyto


def test(name, base_sp, nsp, fs, ttbase, clf=None, clf_feats=None, nfeats=0,
        **kwargs):
    (arr_train, arr_test), npos = ppi.learning_examples(base_sp, fs, ttbase,
            nsp, **kwargs)
    clf_feats = clf_feats if clf_feats else ml.linear()
    if nfeats:
        print 'Selecting top %s of %s features' % (nfeats,
            len(arr_train.dtype.names)-3)
        feats = ml.feature_selection(arr_train, clf_feats)[0][:nfeats]
        arr_train, arr_test = [fe.keep_cols(arr,feats) 
                for arr in arr_train, arr_test]
    else:
        feats = []
    clf = clf if clf else ml.svm()
    scaler = ml.fit_clf(arr_train, clf)
    ppis = ml.classify(clf, arr_test, scaler=scaler)
    result = Struct(train=arr_train[['id1','id2','hit']], clf=clf,
            scaler=scaler, ppis=ppis, npos=npos, name=name, species=base_sp,
            ppi_params=str(clf), feats=feats)
    return result

def cluster(result, fracppis, **kwargs):
    """
    Works on output of either test or predict.
    """
    nppis = int(fracppis*len(result.ppis))
    result.cxs, result.cxppis = cl.filter_c1(result.ppis[:nppis], **kwargs)
    result.cx_params = str(kwargs)
    return result

def predict(name, sp, arr, clf, nsp, nsp_cutoff=0.25):
    ut.savepy(arr, ut.bigd('all_'+name+'.pyd'))
    assert hasattr(clf, 'n_classes_'), "Classifier not yet trained."
    arr_filt, pd_spcounts =  ml.filter_nsp(arr, nsp=nsp,
            cutoff=nsp_cutoff, maybedontfilt=False)
    assert len(arr_filt)>0, "No interactions passed filter."
    #if nsp==1: 
        #assert (arr==arr_filt).all(), 'species filtering issue'
    #else: 
        #ut.savepy(arr_filt, ut.bigd('all_'+name+'_%ssp.pyd' % nsp))
    print "Filtering reduced row count from %s to %s" %(len(arr),len(arr_filt))
    ppis = ml.classify(clf, arr_filt)
    return Struct(ppis=ppis,name=name, species=sp, pd_spcounts=pd_spcounts,
            ppi_params=str(clf))  

def cyto_export(result, arrtrain, ppis=None, cxs=None, geneatts=
        'Hs_ensg_name_desc_uni_entrez.tab', species=None, name_ext=''):
    fname = 'cy_'+result.name+name_ext+'.tab'
    species = species if species else result.species
    ppis = ppis if ppis else result.cxppis
    cxs = cxs if cxs else result.cxs
    cyto.cyto_prep(result.cxppis, arrtrain, fname, geneatts, cxs=result.cxs,
            species=species, pd_spcounts=result.pd_spcounts)

def replace_eluts(arr, scores, name='eluts'):
    species = set(ut.config()['elut_species'].split('_'))
    keep_cols = [n for n in arr.dtype.names[3:] if n[:2] not in species]
    keep_cols = keep_cols + [name]
    newarr = ppi.base_array(arr[['id1','id2','hit']], keep_cols, len(arr))
    newarr[name] = scores
    return newarr
