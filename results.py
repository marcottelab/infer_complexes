import numpy as np
import utils as ut
from Struct import Struct
import myml
import cluster as cl
import ppi
import cyto


def test(name, sp, test_train_npos, clf, nsp, nsp_cutoff=0.25):
    (arr_train, arr_test), npos = test_train_npos
    myml.fit_clf(arr_train, clf)
    arr_test_filt, pd_spcounts, _ =  myml.filter_nsp(arr_test, nsp=nsp,
            cutoff=nsp_cutoff, maybedontfilt=False)
    if nsp==1: assert (arr_test==arr_test_filt).all(), 'species filtering issue'
    ppis = myml.classify(clf, arr_test_filt)
    result = Struct(train=arr_train[['id1','id2','hit']], clf=clf, ppis=ppis, 
            npos=npos, name=name, species=sp, ppi_params=str(clf),
            pd_spcounts=pd_spcounts)
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
    arr_filt, pd_spcounts, _ =  myml.filter_nsp(arr, nsp=nsp,
            cutoff=nsp_cutoff, maybedontfilt=False)
    assert len(arr_filt)>0, "No interactions passed filter."
    #if nsp==1: 
        #assert (arr==arr_filt).all(), 'species filtering issue'
    #else: 
        #ut.savepy(arr_filt, ut.bigd('all_'+name+'_%ssp.pyd' % nsp))
    print "Filtering reduced row count from %s to %s" %(len(arr),len(arr_filt))
    ppis = myml.classify(clf, arr_filt)
    return Struct(ppis=ppis,name=name, species=sp, pd_spcounts=pd_spcounts,
            ppi_params=str(clf))  

def cyto_export(result, arrtrain, ppis=None,
        geneatts='Hs_ensg_name_desc_uni_entrez.tab', species=None):
    fname = 'cy_'+result.name+'.tab'
    species = species if species else result.species
    ppis = ppis if ppis else result.cxppis
    cyto.cyto_prep(result.cxppis, arrtrain, fname, geneatts, cxs=result.cxs,
            species=species, pd_spcounts=result.pd_spcounts)

def double_up(task, name, species, elutfs, clf1, clf2, cutoff1=0.25, nsp=2,
        cutoff2=0.25, **kwargs):
    """
    task: {'test', 'predict'}
    """
    testing = task=='test'
    if testing:
        ppi_func = ppi.learning_examples
        assert 'base_tt' in kwargs, 'Must supply base_tt in test'
    else:
        ppi_func = ppi.predict_all
    test1 = ppi_func(species, elutfs, cutoff=cutoff1, **kwargs)
    if testing:
        (train1, test1), npos = test1
        myml.fit_clf(train1, clf1)
    _, pd_spcounts, orig_spcounts = myml.filter_nsp(test1, nsp=nsp,
            cutoff=cutoff2, maybedontfilt=False)
    def new_feature(clf, to_test):
        # don't sort, since joining these in with array later.
        scores = [p[2] for p in myml.classify(clf, to_test, do_sort=False)]
        return scores
    if testing:
        train2, test2 = [replace_eluts(t, new_feature(clf1, t)) for t in
                train1, test1]
        myml.fit_clf(train2, clf2)
    else:
        test2 = replace_eluts(test1, new_feature(clf1, test1))
    # Filter based on species after doing the fitting--may as well learn from
    # these examples before tossing them. But before classifying to save time.
    test2 = test2[np.array(orig_spcounts) >= nsp]
    ppis_nsp = myml.classify(clf2, test2)
    result = Struct(ppis=ppis_nsp, name=name, clf1=clf1, clf2=clf2,
            ppi_params= str(kwargs)+str(clf1)+str(clf2), pd_spcounts=
            pd_spcounts, ppis_fromfirstclf=test2[['id1','id2','eluts','hit']])
    if testing:
        result.npos = npos
        result.train = kwargs['base_tt'][0]
        result.ppis_fromfirstclf = sorted([p for p in test2[['id1','id2','eluts','hit']]],key=lambda x: x[2], reverse=True)
    return result

def replace_eluts(arr, scores, name='eluts'):
    species = set(ut.config()['elut_species'].split('_'))
    keep_cols = [n for n in arr.dtype.names[3:] if n[:2] not in species]
    keep_cols = keep_cols + [name]
    newarr = ppi.base_array(arr[['id1','id2','hit']], keep_cols, len(arr))
    newarr[name] = scores
    return newarr
