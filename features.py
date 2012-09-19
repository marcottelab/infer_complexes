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

def merge_features(arr, pattern, func, do_remove):
    feats = ut.regex_filter(arr.dtype.names, pattern)
    print "features matched:", feats
    merged = [func(i) for i in arr[feats]]
    name = pattern + '_' + func.__name__
    removefeats = set(feats) if do_remove else set([])
    keepdtypes = [dt for dt in arr.dtype.descr if not dt[0] in removefeats]
    newdtype = np.dtype(keepdtypes + [(name, arr[feats[0]].dtype)])
    newarr = np.empty(arr.shape, dtype=newdtype)
    for field in [d[0] for d in keepdtypes]:
        newarr[field] = arr[field]
    newarr[name] = merged
    return newarr

def retype_arr(arr, oldtype='f2', newtype='f4'):
    newdtype = [(dt[0],dt[1].replace(oldtype,newtype)) 
            for dt in arr.dtype.descr]
    newarr = np.empty(arr.shape, dtype=newdtype)
    for field in newarr.dtype.names:
        newarr[field] = arr[field]
    return newarr

def merge_by_species(arr, matches, func, remove=False):
    """
    matches: like [apex] or [wcc, apex, ...]
    Makes patterns with match for each species, like 'Hs.*apex
    """
    assert not isinstance(matches,str), "matches is list, not string"
    def merge_recurse(arr, patterns, func):
        if patterns:
            newarr = merge_features(arr, patterns[0], func, remove)
            return merge_recurse(newarr, patterns[1:], func)
        else:
            return arr
    # Won't match a merged feature since that will have a * in it.
    patterns = [sp+'\w*'+match for match in matches for sp in
            ut.config()['elut_species'].split('_')]
    return merge_recurse(arr, patterns, func)


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

def filter_require_sp(arr, species, cutoff=0.25, count_ext=True):
    features = arr.dtype.names[3:]
    cols = [f for f in features if f[:2]==species or (count_ext==True and
        f[:2]=='ex' and f[4:6]==species) ]
    sp_max = [max(r) for r in arr[cols]]
    passing_inds = [i for i,m in enumerate(sp_max) if m > cutoff]
    return arr[passing_inds]

def filter_nsp(arr, nsp=2, cutoff=0.25, maybedontfilt=True, count_ext=True):
    """
    Filter to only leave interactions for which evidence exists in nsp species.
    scores: [arr_train, arr_test]
    In addition to returning the filtered array, returns a matching pairdict of
    species counts.
    """
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
    sp_counts_filt = [c for c in orig_sp_counts if c >=nsp]
    arrfilt = arr[exceed_inds]
    pd_spcounts = pd.PairDict([[arrfilt[i][0],arrfilt[i][1],sp_counts_filt[i]] 
                            for i in range(len(arrfilt))])
    return arrfilt, pd_spcounts

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
