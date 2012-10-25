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
from sklearn.svm import LinearSVC
import utils as ut
import numpy as np
import multiprocessing
import features as fe

NCORES = multiprocessing.cpu_count()

def fit_and_test(scored, clf, norm=True): 
    """
    scored: [arr_train, arr_test]
    """
    arr_train, arr_test = scored
    if not exist_pos_neg(arr_train):
        return []
    scaler = fit_clf(arr_train, clf, norm=norm)
    tested = classify(clf, arr_test, scaler)
    return tested

def exist_pos_neg(arr):
    hits = arr['hit']
    found_true = False
    found_false = False
    for h in hits:
        if h:
            found_true = True
        else:
            found_false = True
        if found_true and found_false:
            return True
    return False

def fit_clf(arr, clfbase, norm=True):
    if norm:
        arr = fe.retype_arr(arr) # change f2 to f4 to prevent overflow
    X,y = arr_feats(arr), arr['hit']
    scaler = None
    if norm:
        print "Fitting and scaling training features."
        scaler = sk.preprocessing.Scaler().fit(X)
        X = scaler.transform(X)
    print "Training classifier: %s examples, %s features" % (len(X), len(X[0]))
    clfbase.fit(X,y)
    return scaler

def classify(clf, arr, scaler=None, do_sort=True):
    """
    If the clf was trained without at least a pos and a neg, this will fail.
    """
    X = arr_feats(arr)
    if scaler: 
        print "Scaling features before prediction."
        X = scaler.transform(X)
    print "Predicting: %s examples, %s features" % (len(arr), len(X[0]))
    probs = (x[1] for x in clf.predict_proba(X))
    tested = zip(arr['id1'], arr['id2'], probs, arr['hit'])
    if do_sort: tested.sort(key=lambda x:x[2],reverse=True)
    return tested

def arr_feats(arr):
    return [[x for x in r] for r in arr[list(arr.dtype.names[3:])]]

def tree(n_estimators=200,n_jobs=NCORES-1, bootstrap=True, **kwargs):
    return ExtraTreesClassifier(n_estimators=n_estimators, n_jobs=n_jobs,
                                bootstrap=bootstrap, **kwargs)

def svm(kernel='linear', cache_size=4000, **kwargs):
    return SVC(kernel=kernel, cache_size=cache_size, probability=True, **kwargs)

def linear(dual=False, **kwargs):
    return LinearSVC(dual=dual, **kwargs)

def feature_selection(arr, clf, do_plot=False):
    """
    clf: ml.tree(compute_importances=True) or ml.linear()
    """
    names = arr.dtype.names[3:]
    fit_clf(arr, clf, norm=True)
    importances = (clf.coef_[0] if hasattr(clf, 'coef_') else
            clf.feature_importances_)
    indices = np.argsort(importances)[::-1]
    ranked = [(names[index], importances[index]) for index in indices]
    for i,(name,imp) in enumerate(ranked):
        print "%d. %s (%f)" % (i + 1, name, imp)
    # Plot the feature importances of the trees and of the forest
    if do_plot:
        import pylab as pl
        pl.figure()
        pl.title("Feature importances")
        for tree in forest.estimators_:
            pl.plot(indnums, tree.feature_importances_[indices], "r")
        pl.plot(indnums, importances[indices], "b")
        pl.show()
    feats, weights = zip(*ranked)
    return list(feats), list(weights)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("usage: python ml.py train_test feats_f clf_type \
               donorm kwarg1_val1-kwarg2-val2")
    ttf = sys.argv[1]
    tt = np.load(ttf)
    feats = ut.loadpy(sys.argv[2])
    k = sys.argv[3]
    do_norm = sys.argv[4]
    kvs = sys.argv[5]
    kwargs = dict([tuple(kv.split('_')) for kv in kvs.split('-')]) \
        if kvs else {}
    clf = tree(**kwargs) if k=='tree' else svm(kernel=k, **kwargs)
    ts =  [('%s features, %s kernel, norm: %s, %s' %(n,k,do_norm, kvs),
        fit_and_test([fe.keep_cols(t, ut.i0(feats[:n])) for t in tt], 
                        clf, norm=do_norm)) 
        for n in 20,30,40,50]
    ut.savepy(ts, 'ts_%s_%s_%s_%s' %(k,do_norm,kvs,ttf))
