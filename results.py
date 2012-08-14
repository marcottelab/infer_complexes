import utils as ut
from Struct import Struct
import myml
import cluster as cl
import ppi

def test(name, test_train_npos, clf, **kwargs):
    (arr_train, arr_test), npos = test_train_npos
    myml.fit_clf(arr_train, clf, **kwargs)
    ppis = myml.classify(clf, arr_test, **kwargs)
    result = Struct(train=arr_train, test=arr_test, clf=clf, ppis=ppis, 
            npos=npos, name=name, ppi_params=str(kwargs)+str(clf))
    return result

def cluster(result, fracppis, **kwargs):
    """
    Works on output of either test or predict.
    """
    nppis = int(fracppis*len(result.ppis))
    result.cxs, result.cxppis = cl.filter_c1(result.ppis[:nppis], **kwargs)
    result.cx_params = str(kwargs)
    return result

def predict(name, arr, clf, **kwargs):
    ppis = myml.classify(clf, arr, **kwargs)
    return Struct(ppis=ppis,name=name)

"""
rename to pipeline.  take over the running of ppi.  run it the first time as
is, and the second call to either all_predictions or learning_examples removing
elut_fs and providing the output scores as other_data.
"""
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
    print '>>test1 len', len(test1)
    trash_arr, pd_spcounts = myml.filter_nsp(test1, nsp=nsp, cutoff=cutoff2,
            dontfilt=False)
    print '>>spcounts len', len(pd_spcounts.d)
    new_feature = (p[:3] for p in myml.classify(clf1, test1))
    
    test2 = ppi_func(species, [], cutoff=-1,
            other_evidence=[('eluts', new_feature)], **kwargs)
    if testing:
        (train2, test2), npos = test2
        myml.fit_clf(train2, clf2)
    print '>>test2 len', len(test2)
    # Filter based on species after doing the fitting--may as well learn from
    # these examples before tossing them. But before classifying to save time.
    test2 = filter_nsp_pd(test2, nsp, pd_spcounts)
    print '>>test2 filtered len', len(test2)
    ppis_nsp = myml.classify(clf2, test2)
    result = Struct(ppis=ppis_nsp, name=name, clf1=clf1, clf2=clf2,
            ppi_params= str(kwargs)+str(clf1)+str(clf2), pd_spcounts=
            pd_spcounts)
    if testing:
        result.npos = npos
        result.trainexs = kwargs[base_tt[0]]
    return result
