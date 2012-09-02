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
from ppi import filter_arr
import numpy as np
import multiprocessing
import operator
import pairdict as pd

NCORES = multiprocessing.cpu_count()
NET_SPS = 'HS CE DM SC'.split()

def boot_arr(arr):
    return arr[ut.sample_wr(range(len(arr)), len(arr))]

def match_cols(arr, colstring):
    """
    colstring like 'Hs Ce HS ...' or 'Hs Ce HS-GN ...'
    """
    return [f for f in arr.dtype.names if (f[:2] in colstring.split()
                                           or f in colstring.split())]
def keep_cols(arr, colnames):
    return arr[['id1','id2','hit'] + colnames]

def regex_cols(arr, pattern):
    return arr[['id1','id2','hit'] + ut.regex_filter(arr.dtype.names, pattern)]

def merge_features(arr, pattern, func):
    feats = ut.regex_filter(arr.dtype.names, pattern)
    merged = reduce(func, [arr[feat] for feat in feats])
    name = pattern + '_' + func.__name__
    setfeats = set(feats)
    keepdtypes = [dt for dt in arr.dtype.descr if not dt[0] in setfeats]
    newdtype = np.dtype(keepdtypes + [(name, merged.dtype)])
    newarr = np.empty(arr.shape, dtype=newdtype)
    for field in [d[0] for d in keepdtypes]:
        newarr[field] = arr[field]
    newarr[name] = merged
    return newarr

def merge_by_species(arr, match, func):
    """
    match: like wcc, apex, ...
    Makes patterns with match for each species, like 'Hs.*apex
    """
    def merge_recurse(arr, patterns, func):
        if patterns:
            newarr = merge_features(arr, patterns[0], func)
            return merge_recurse(newarr, patterns[1:], func)
        else:
            return arr
    patterns = [sp+'.*'+match for sp in ut.config()['elut_species'].split('_')]
    return merge_recurse(arr, patterns, func)

def feature_selection(arr, columns=None, cutoff=0.25, n_est=50,
        n_jobs=NCORES-1, do_plot=False):
    # Build a forest and compute the feature importances
    forest = ExtraTreesClassifier(n_estimators=n_est,
                                  compute_importances=True,
                                  random_state=0, n_jobs=n_jobs)
    arr,names = filter_names_arr(arr, columns, cutoff)
    X = features(arr, names)
    y = arr['hit']
    forest.fit(X, y)
    importances = forest.feature_importances_
    indices = np.argsort(importances)[::-1]
    # Print the feature ranking
    ranked = [(names[index], importances[index]) for index in indices]
    print "Feature ranking:"
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
    return ranked

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

def filter_nsp(arr, nsp=2, cutoff=0.25, maybedontfilt=True, count_ext=True):
    """
    Filter to only leave interactions for which evidence exists in nsp species.
    scores: [arr_train, arr_test]
    In addition to returning the filtered array, returns a matching pairdict of
    species counts.
    """
    if maybedontfilt and nsp==1 and cutoff==0.25:
        print "Assuming no need to filter: nsp=%s, cutoff=%s" % (nsp, cutoff)
        return arr, None
    if len(arr) == 2:
        # user provided [trainarr, testarr]
        return [filter_nsp(a, nsp=nsp, cutoff=cutoff,
            maybedontfilt=maybedontfilt, count_ext=count_ext)
                for a in arr]
    features = arr.dtype.names[3:]
    assert 'eluts' not in features, "Counting sps not supported with 'eluts'"
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
    orig_sp_counts = [len([1 for m in maxes if m[i] > cutoff]) 
        for i in range(len(arr))]
    sp_counts_filt = [c for c in orig_sp_counts if c >=nsp]
    arrfilt = arr[exceed_inds]
    pd_spcounts = pd.PairDict([[arrfilt[i][0],arrfilt[i][1],sp_counts_filt[i]] 
                            for i in range(len(arrfilt))])
    return arrfilt, pd_spcounts, orig_sp_counts

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
