import utils as ut
from Struct import Struct
import myml
import cluster as cl

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
