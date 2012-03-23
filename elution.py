from __future__ import division
import numpy as np
import operator
import random
import os
from Struct import Struct
import utils as ut
import corr
import cv

def load_elution_desc(fname):
    # expected file structure:
    # first col: gene id
    # second col: gene description
    # remaining cols: elution profile data
    mat = np.matrix([row[2:] for row in ut.load_tab_file(fname)][1:],dtype='float64')
    (prots,gdesc) = zip(*[(row[0],row[1]) for row in ut.load_tab_file(fname)][1:])
    elut = Struct(mat=mat, prots=prots, gdesc=gdesc, filename=fname)
    return elut

def load_elution_total(fname):
    # expected file structure:
    # first col: gene id
    # second col: total
    # remaining cols: elution profile data
    # final row: total count
    rows = [r for r in ut.load_tab_file(fname)][:-1]
    mat = np.matrix([row[2:] for row in rows][1:],dtype='float64')
    (prots,totals) = zip(*[(row[0],row[1]) for row in rows][1:])
    elut = Struct(mat=mat, prots=prots, totals=totals, filename=fname,
                  filename_original=fname)
    return elut

def load_complexes(filename, format_single_protein=False):
    # load corum-type file into a dictionary
    # complexes: dict{complexid: set([protein1, protein2,...]), .. }
    # first col: complex id
    # third col: protein id
    if format_single_protein:
        complexes = {}
        for l in ut.load_tab_file(filename):
            complexes.setdefault(l[0],set([])).add(l[2])
    else:
        complexes = dict([(l[0],set(l[1:])) for l in
                          ut.load_list_of_lists(filename)])
    return complexes
            
    
def load_interactions(filename):
    # complexes: dict{complexid: set([protein1, protein2,...]), .. }
    # this makes a dictionary{protein1: set([protein2, protein3]), ...}
    # every interaction is found twice here for fast interaction checking
    complexes = load_complexes(filename)
    interactions = {}
    for complex,protein_set in complexes.items():
        for p in protein_set:
            partners = protein_set.copy()
            partners.remove(p)
            [interactions.setdefault(p,set([])).add(par) for par in partners]
    return interactions

def process_raw_wan(f_source, f_dest=None, first_col_element=1,
                    first_data_col=1, end_description_col=True,
                    first_data_row=1):
    # specific to cuihong's files, and tries to handle the differences seen in them
    # always keeps the first column as variable name
    # processes first column, splitting it and keeping first_col_element
    # for the array, keeps columns [first_data_col:end_data_col]. None works.
    # textlines = [textline for textline in open(f_source)]
    # # handle unix-unreadable linebreaks from excel
    # if len(textlines) == 1:
    #     if textlines[0].find('\r\n') > -1:
    #         textlines = textlines[0].split('\r\n')
    lines = [line.strip().split('\t') for line in open(f_source)if line.strip()!='']
    # simple: one step at a time.
    # column manipulation first.
    if end_description_col:
        lines = [[l[0]] + [l[-1]] + l[first_data_col:-1] for l in lines]
    else:
        lines = [[l[0]] + l[first_data_col:] for l in lines]
    # variable name manipulation
    if first_col_element is not None:
        # manipulate gene name in all but header row. skip anything btw header
        # and first_data_row.
        lines = [lines[0]] + [[l[0].split('|')[first_col_element]] +
                    l[1:] for l in lines[first_data_row:]]
    # rename file
    if f_dest is None:
        split = os.path.splitext(f_source)
        f_dest = split[0] + '_proc' + split[1]
    ut.write_tab_file(lines, f_dest)

def correlate_single(elut1, elut2, prot):
    return np.corrcoef(elut1.mat[elut1.prots.index(prot),:],
        elut2.mat[elut2.prots.index(prot),:])[0][1]

def correlate_matches(elut1, elut2):
    overlap = set.intersection(set(elut1.prots),set(elut2.prots))
    return [correlate_single(elut1, elut2, p) for p in overlap]
    
def correlate_matches_dict(elut1, elut2, pdict_1to2):
    overlap = [p for p in elut1.prots if p in pdict_1to2 and
        list(pdict_1to2[p])[0] in set(elut2.prots)]
    return [(p,np.corrcoef(elut1.mat[elut1.prots.index(p),:],
        elut2.mat[elut2.prots.index(list(pdict_1to2[p])[0]),:])[0][1]) for p in
        overlap]
    
def compute_cvpairs(elutfile, trueints, sample_frac, poisson_repeat=200):
    basename = os.path.splitext(os.path.split(elutfile)[1])[0]
    elut = load_elution_total(elutfile)
    elut.corr = corr.traver_corr(elut.mat, repeat=poisson_repeat)
    return cv.cv_pairs(elut.corr, trueints, elut.prots, sample_frac=sample_frac)

def combine_elutions(e1, e2, combine_corr_func=None):
    # functions: np.maximum, sum, ...
    allprots = list(set.union(set(e1.prots), set(e2.prots)))
    nprots = len(allprots)
    allfracs = e1.mat.shape[1] + e2.mat.shape[1]
    mat = np.matrix(np.zeros((nprots,allfracs)))
    mat[0:len(e1.prots),0:e1.mat.shape[1]] = e1.mat[:,:]
    for row2 in range(len(e2.prots)):
        mat[allprots.index(e2.prots[row2]), e1.mat.shape[1]:] = e2.mat[row2,:]
    elut = Struct(mat=mat, prots=allprots,
                  filename=e1.filename+e2.filename+str(combine_corr_func))
    if combine_corr_func:
        elut.corr = combine_corrs(e1, e2, allprots, combine_corr_func)
    return elut

def combine_corrs(e1, e2, allprots, combine_func, default_val=None):
    # we combine the symmetric correlation matrices using the specified
    # element-wise function. function examples: max, sum
    # we use the specified ordering of elements in allprots
    default_val = default_val if default_val else -1 if \
        combine_func.__name__.find('max') > -1 else 0
    nprots = len(allprots)
    corr = np.matrix(np.zeros((nprots,nprots)))
    dprots1 = ut.list_inv_to_dict(e1.prots)
    dprots2 = ut.list_inv_to_dict(e2.prots)
    for row,p1 in enumerate(allprots):
        for col,p2 in enumerate(allprots):
            val1 = e1.corr[dprots1[p1], dprots1[p2]] if p1 in dprots1 and p2 in \
                dprots1 else default_val
            val2 = e1.corr[dprots2[p1], dprots2[p2]] if p1 in dprots2 and p2 in \
                dprots2 else default_val
            corr[row,col] = combine_func(val1, val2)
    return corr
    
def test_combined_corrs(eluts, ncomparisons=10):
    # compare at ncomparisons randomly-selected places
    prots_common = list(reduce(set.union,[set(e.prots) for e in eluts]))
    p1s = random.sample(prots_common, ncomparisons)
    p2s = random.sample(prots_common, ncomparisons)
    return [[e.corr[e.prots.index(p1),e.prots.index(p2)] for e in eluts] for
(p1,p2) in zip(p1s, p2s)]
