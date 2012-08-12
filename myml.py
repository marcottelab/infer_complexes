from __future__ import division
import os
import sys
if not '/home/blakeb/.local/lib/python2.7/scikit_learn-0.11-py2.7-linux-x86_64.egg' in sys.path:
    if os.path.exists('/home/blakeb/.local/lib/python2.7/scikit_learn-0.11-py2.7-linux-x86_64.egg'):
        sys.path.append('/home/blakeb/.local/lib/python2.7/scikit_learn-0.11-py2.7-linux-x86_64.egg')
from sklearn.svm import SVC
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
import utils as ut
from ppi import filter_arr
import numpy as np
import multiprocessing
import operator
import pairdict as pd

NCORES = multiprocessing.cpu_count()
NET_SPS = 'HS CE DM SC'.split()

def fit_and_test(scored, clf, columns=None, cutoff=0.25): 
    """
    scored: [arr_train, arr_test]
    """
    arr_train, arr_test = scored
    fit_clf(arr_train, clf, columns)
    tested = classify(clf, arr_test, columns)
    return tested

def fit_clf(arr_train, clfbase, columns=None, cutoff=0.25):
    arr_train, names = filter_names_arr(arr_train, columns, cutoff)
    print "Training classifier: %s examples, %s features" % (len(arr_train),
        len(names))
    clfbase.fit(features(arr_train, names), arr_train['hit'])

def classify(clf, arr_test, columns=None, cutoff=0.25, savef=None):
    arr_test, names = filter_names_arr(arr_test, columns, cutoff)
    if columns: print 'Features:', ut.count_collect(names, 6)
    print "Classifying: %s examples, %s features" % (len(arr_test), 
            len(names))
    probs = (x[1] for x in clf.predict_proba(features(arr_test, names)))
    tested = zip(arr_test['id1'], arr_test['id2'], probs, arr_test['hit'])
    tested.sort(key=lambda x:x[2],reverse=True)
    if savef: ut.savepy(tested, savef)
    return tested

def classify_slice(clf,arr_test, perslice, savef=None, columns=None,
        cutoff=0.25, maintain=True, startslice=0):
    nslices = int(np.ceil(len(arr_test) / perslice))
    slices = (classify(clf, arr_test[i*perslice:(i+1)*perslice],
                       savef=(savef+str(i)+'.pyd') if savef else None,
                       columns=columns, cutoff=cutoff)
              for i in range(startslice, nslices))
    if maintain==True:
        tested = reduce(operator.add, slices)
        tested.sort(key=lambda x:x[2],reverse=True)
        return tested
    else:
        for i,s in enumerate(slices):
            print i, savef

def tree(n_estimators=100,n_jobs=int(NCORES/2), bootstrap=True, **kwargs):
    return ExtraTreesClassifier(n_estimators=n_estimators, n_jobs=n_jobs,
                                bootstrap=bootstrap, **kwargs)

def svm(kernel='linear', prob=True, **kwargs):
    return SVC(kernel=kernel, probability=prob, **kwargs)

def filter_names_arr(arr, columns, cutoff, nofilter=set(NET_SPS)):
    """
    columns: either a list of column numbers or a space-sep string of
    2-letter matches for column names.
    Filters the given array, returning the full rows (not just those named
    columns) for which at least one item passes the given cutoff.
    Those beginning with items in nofilter are not used in qualifying
    threshold-passing rows.
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
        all_possible = [n for n in arr.dtype.names if not n[:2] in nofilter]
        if set(names_filt) == set(all_possible):
            print 'no filtering'
            newarr = arr
        else:
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

def filter_nsp(arr, nsp=2, cutoff=0.25, dontfilt=True, count_ext=True):
    """
    Filter to only leave interactions for which evidence exists in nsp species.
    scores: [arr_train, arr_test]
    In addition to returning the filtered array, returns a matching pairdict of
    species counts.
    """
    if len(arr) == 2:
        # user provided [trainarr, testarr]
        return [filter_nsp(a, nsp=nsp, cutoff=cutoff, dontfilt=dontfilt,
            count_ext=count_ext)
                for a in arr]
    if dontfilt and nsp==1 and cutoff==0.25:
        print "Assuming no need to filter: nsp=%s, cutoff=%s" % (nsp, cutoff)
        return arr
    features = arr.dtype.names[3:]
    sps = set([n[:2] for n in features 
                if n[:2] not in set(NET_SPS) and n[:2] != 'ex'])
    print 'Filtering >=%s of these species >=%s:' % (nsp, cutoff), sps
    spcols = [[f for f in features 
                if f[:2]==s or (count_ext==True and f[:2]=='ex' and f[4:6]==s)]
                for s in sps]
    # Note that maxes are transposed relative to what might be expected.
    maxes = [[max(i) for i in arr[scs]] for scs in spcols]
    exceed_inds = [i for i in range(len(arr))
                   if len([i for m in maxes if m[i] > cutoff]) >= nsp]
    sp_counts = [c for c in [len([1 for m in maxes if m[i] > cutoff]) 
        for i in range(len(arr))] if c >= nsp]
    arrfilt = arr[exceed_inds]
    pd_spcounts = pd.PairDict([[arrfilt[i][0],arrfilt[i][1],sp_counts[i]] 
                            for i in range(len(arrfilt))])
    return arr[exceed_inds], pd_spcounts

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 6:
        sys.exit("usage: python myml.py f_examples f_classifier \
               perslice islice path")
    exs = np.load(sys.argv[1])
    clf = ut.loadpy(sys.argv[2])
    perslice = int(sys.argv[3])
    i = int(sys.argv[4])
    path = sys.argv[5]
    exs_slice = exs[i*perslice:(i+1)*perslice]
    del exs
    classify_slice(clf, exs_slice, 100000, savef=(path+str(i)+'_'),
            maintain=False, startslice=0)
