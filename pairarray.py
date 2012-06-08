from __future__ import division
import utils as ut
from numpy import ndarray
import pairdict as pd

class PairArray(object):
    """
    Has a PairDict that indexes to rows of the array.
    self.names: Stores a name for each column in the array.
    Arr[:,0] is {1,0,-1} for t,f,unknown interaction: 
    """
    
    def __init__(self, pdict, arr, names):
        """
        Does NOT copy arr.
        Copies names, pd.
        """
        self.array = arr
        self.names = list(names)
        self.pdict = pd.pd_copy(pdict)
