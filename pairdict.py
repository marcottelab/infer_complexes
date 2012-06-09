from __future__ import division
import utils as ut
from numpy import ndarray

class PairDict(object):
    """
    A common storage format for examples and predictions that handles all the
    messiness of merging, deduping, etc.
    """

    def __init__(self, lopair_vals):
        """
        Make this from a list or set of tuples of (id1, id2, val1, val2, ...)
        """
        self.d = dict([((p[0],p[1]),list(p[2:])) for p in lopair_vals])

    def set(self, key, val):
        k = self.find(key)
        if k==None:
            k = key
        self.d[k] = val

    def find(self, pair):
        if pair in self.d:
            return pair
        elif pd_flip(pair) in self.d:
            return pd_flip(pair)
        else:
            return None

def pd_copy(pd):
    newpd = PairDict([])
    newpd.d = pd.d.copy()
    return newpd

def pd_flip(pair):
    return (pair[1],pair[0])

def pd_lol(pd):
    return [[k[0],k[1]] + pd.d[k] for k in pd.d]

def pd_union(a,b,adefaults=None,bdefaults=None):
    """
    Merge two PairDicts and return a new one.
    Values for each pair become a list of avalues+bvalues, with defaults for
    each index provided if desired.
    """
    adefaults = adefaults if adefaults else [None]*len(a.d.values()[0])
    bdefaults = bdefaults if bdefaults else [None]*len(b.d.values()[0])
    newpd = PairDict([])
    def merge_help(from_set,newpd,a,b,adefaults,bdefaults, reverse=False):
        bleftovers = set(b.d.keys())
        for apair in from_set:
            bpair = b.find(apair)
            if bpair:
                newvals = (a.d[apair] + b.d[bpair] if not reverse else
                            b.d[bpair] + a.d[apair])
                bleftovers.remove(bpair)
            else:
                newvals = (a.d[apair] + bdefaults if not reverse else bdefaults
                        + a.d[apair])
            newpd.set(apair,newvals)
        return bleftovers
    bleftovers = merge_help(a.d.keys(),newpd,a,b,adefaults,bdefaults)
    merge_help(bleftovers,newpd,b,a,bdefaults,adefaults,reverse=True)
    return newpd
