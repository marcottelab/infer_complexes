from sklearn import svm
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
import utils as ut
from ppi import filter_arr
import numpy as np

def fit_and_test(scored, columns=None, cutoff=0.25): #[arr_train, arr_test]
    arr_train, arr_test = scored
    classer = fit_svm(arr_train, columns)
    tested = test_svm(classer, arr_test, columns)
    return tested

def fit_svm(arr_train, columns=None, cutoff=0.25):
    classer = svm.SVC(kernel='linear', probability=True) 
    arr_train, names = filter_names_arr(arr_train, columns, cutoff)
    classer.fit(features(arr_train, names), arr_train['hit'])
    return classer

def test_svm(classer, arr_test, columns=None, cutoff=0.25):
    arr_test, names = filter_names_arr(arr_test, columns, cutoff)
    print 'Features:', names
    probs = (x[1] for x in classer.predict_proba(features(arr_test, names)))
    tested = zip(arr_test['id1'], arr_test['id2'], probs, arr_test['hit'])
    tested.sort(key=lambda x:x[2],reverse=True)
    return tested

def filter_names_arr(arr, columns, cutoff, nofilter=['HS','CE','DM','SC']):
    """
    columns: either a list of column numbers or a space-sep string of
    2-letter matches for column names.
    Filters the given array, returning the full rows (not just those named
    columns) for which at least one item passes the given cutoff.
    """
    columns = columns if columns else feature_inds(arr)
    feat_names = ([arr.dtype.names[i] for i in columns]
                  if isinstance(columns, list)
                  else match_cols(arr,columns))
    feat_nums = ([i for i,name in enumerate(arr.dtype.names)
                  if name in feat_names])
    if feat_nums == feature_inds(arr):
        print 'no filtering'
        newarr = arr
    else:
        # DON'T filter by network score columns: these don't qualify the row.
        names_filt = [n for n in feat_names if not n[:2] in (nofilter)]
        print 'filtering', names_filt
        nums_filt = ([i for i,name in enumerate(arr.dtype.names)
                  if name in names_filt])
        newarr = filter_arr(arr, nums_filt, cutoff)
    return newarr, feat_names

def features(arr, feat_names):
    return [[x for x in l] for l in arr[feat_names]]

def feature_inds(arr):
    return range(3, len(arr[0]))

def match_cols(arr, colstring):
    """
    colstring like 'Hs Ce HS ...' or 'Hs Ce HS-GN ...'
    """
    return [f for f in arr.dtype.names if (f[:2] in colstring.split()
                                           or f in colstring.split())]

def rescale(p,n):
    """
    Rescale posterior probability p according to multiple of negatives n.
    """
    return p/(1+(1-p)*n)

def feature_selection(arr, columns=None, cutoff=0.25):

    # Build a forest and compute the feature importances
    forest = ExtraTreesClassifier(n_estimators=250,
                                  compute_importances=True,
                                  random_state=0)
    arr,names = filter_names_arr(arr, columns, cutoff)
    X = features(arr, names)
    y = arr['hit']
    forest.fit(X, y)
    importances = forest.feature_importances_
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print "Feature ranking:"

    indnums = xrange(len(indices))
    for f in indnums:
        print "%d. %s (%f)" % (f + 1, names[indices[f]],
    importances[indices[f]])
        

    # Plot the feature importances of the trees and of the forest
    import pylab as pl
    pl.figure()
    pl.title("Feature importances")
    for tree in forest.estimators_:
        pl.plot(indnums, tree.feature_importances_[indices], "r")

    pl.plot(indnums, importances[indices], "b")
    pl.show()
