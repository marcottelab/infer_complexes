from __future__ import division
import os
import sys
if not '/home/blakeb/.local/lib/python2.7/scikit_learn-0.11-py2.7-linux-x86_64.egg' in sys.path:
    if os.path.exists('/home/blakeb/.local/lib/python2.7/scikit_learn-0.11-py2.7-linux-x86_64.egg'):
        sys.path.append('/home/blakeb/.local/lib/python2.7/scikit_learn-0.11-py2.7-linux-x86_64.egg')
import sklearn as sk
from sklearn.svm import SVC
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
import utils as ut
import numpy as np
import multiprocessing

NCORES = multiprocessing.cpu_count()

def fit_and_test(scored, clf, norm=True): 
    """
    scored: [arr_train, arr_test]
    """
    arr_train, arr_test = scored
    scaler = None
    if norm:
        scaler = sk.preprocessing.Scaler()
        arr_train = 

    fit_clf(arr_train, clf)
    tested = classify(clf, arr_test)
    return tested

def fit_clf(arr, clfbase):
    feats = feature_names(arr)
    print "Training classifier: %s examples, %s features" % (len(arr),
            len(feats))
    X = arr2lol(arr[feats])
    y = arr['hit']
    clfbase.fit(X,y)

def arr2lol(arr):
    return [[x for x in r] for r in arr]

def feature_names(arr):
    return arr.dtype.names[3:]

def classify(clf, arr, do_sort=True, feats=None):
    feats = feature_names(arr)
    print "Classifying: %s examples, %s features" % (len(arr), 
            len(feats))
    X = arr2lol(arr[feats])
    probs = (x[1] for x in clf.predict_proba(X))
    tested = zip(arr['id1'], arr['id2'], probs, arr['hit'])
    if do_sort: tested.sort(key=lambda x:x[2],reverse=True)
    return tested

def tree(n_estimators=200,n_jobs=NCORES-1, bootstrap=True, **kwargs):
    return ExtraTreesClassifier(n_estimators=n_estimators, n_jobs=n_jobs,
                                bootstrap=bootstrap, **kwargs)

def svm(kernel='linear', cache_size=2000, **kwargs):
    return SVC(kernel=kernel, probability=True, **kwargs)

