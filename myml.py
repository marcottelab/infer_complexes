from sklearn import svm
import utils as ut

def fit_and_test(scored): #[arr_train, arr_test]
    arr_train, arr_test = scored
    classer = fit_svm(arr_train)
    tested = test_svm(classer, arr_test)
    return tested

def fit_svm(arr_train):
    classer = svm.SVC(kernel='linear', probability=True) 
    classer.fit(features(arr_train), arr_train['hit'])
    return classer

def test_svm(classer, arr_test):
    probs = (x[1] for x in classer.predict_proba(features(arr_test)))
    tested = zip(arr_test['id1'], arr_test['id2'], probs, arr_test['hit'])
    tested.sort(key=lambda x:x[2],reverse=True)
    return tested

def features(arr):
    features = list(arr.dtype.names[3:])
    lol = [[x for x in l] for l in arr[features]]
    return lol

def rescale(p,n):
    """
    Rescale posterior probability p according to multiple of negatives n.
    """
    return p/(1+(1-p)*n)
