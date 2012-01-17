from Struct import Struct
import utils as ut
import numpy as np

def load_elution(fname='data/elution.tab'):
    arr = np.array([row[3:] for row in ut.load_tab_file(fname)],dtype='float64')
    genes = [row[0] for row in ut.load_tab_file(fname)]
    elut = Struct(arr=arr, genes=genes)
    return elut
