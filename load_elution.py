import numpy as np
import os
from Struct import Struct
import utils as ut


def load_elution(fname='data/elution.tab'):
    # expected file structure:
    # first col: gene id
    # second col: gene description
    # remaining cols: elution profile data
    arr = np.array([row[3:] for row in ut.load_tab_file(fname)][1:],dtype='float64')
    (genes,gdesc) = zip(*[(row[0],row[1]) for row in ut.load_tab_file(fname)][1:])
    elut = Struct(arr=arr, genes=genes, gdesc=gdesc)
    return elut

def load_multi(elut_files):
    # elut_files: [('Hs', 'folder/filename.tab'), ..]
    eluts = dict([(name,load_elution(fname)) for name,fname in elut_files])
    names = [x[0] for x in elut_files]
    multi_elut = Struct(names=names, eluts=eluts)
    return multi_elut

def process_raw_wan(f_source, f_dest=None, first_col_element=1,
                    first_data_col=1, end_description_col=True):
    # specific to cuihong's files, and tries to handle the differences seen in them
    # always keeps the first column as variable name
    # processes first column, splitting it and keeping first_col_element
    # for the array, keeps columns [first_data_col:end_data_col]. None works.
    lines = [line.strip().split('\t') for line in open(f_source) if line.strip()!='']
    # simple: one step at a time.
    # column manipulation first.
    if end_description_col:
        lines = [[l[0]] + [l[-1]] + l[first_data_col:-1] for l in lines]
    else:
        lines = [[l[0]] + l[first_data_col:] for l in lines]
    # variable name manipulation
    if first_col_element is not None:
        # skip header row.
        lines = [lines[0]] + [[l[0].split('|')[first_col_element]] + l[1:] for
                    l in lines[1:]]
    # rename file
    if f_dest is None:
        split = os.path.splitext(f_source)
        f_dest = split[0] + '_proc' + split[1]
    ut.write_tab_file(lines, f_dest)
