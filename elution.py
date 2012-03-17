from __future__ import division
import numpy as np
import operator
import os
from Struct import Struct
import utils as ut


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
    
def load_elution_desc(fname):
    # expected file structure:
    # first col: gene id
    # second col: gene description
    # remaining cols: elution profile data
    mat = np.matrix([row[2:] for row in ut.load_tab_file(fname)][1:],dtype='float64')
    (prots,gdesc) = zip(*[(row[0],row[1]) for row in ut.load_tab_file(fname)][1:])
    elut = Struct(mat=mat, prots=prots, gdesc=gdesc)
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
    elut = Struct(mat=mat, prots=prots, totals=totals)
    return elut

def load_multi(elut_files):
    # elut_files: [('Hs', 'folder/filename.tab'), ..]
    eluts = dict([(name,load_elution(fname)) for name,fname in elut_files])
    names = [x[0] for x in elut_files]
    multi_elut = Struct(names=names, eluts=eluts)
    return multi_elut

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

def coelution(elut, norm=1):
    # elut.mat: float matrix of spectral counts
    # elut.genes: gene labels for the rows of the array
    if norm:
        # stupid simple: pearson correlation matrix
        corr = np.corrcoef(elut.mat) # between -1 and 1
    else:
        corr = np.cov(elut.mat)
    return corr

def traver_corr(mat, repeat=1000, norm=True):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each row to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_corr(mat, norm=True):
        M = mat.shape[1]
        C = mat + 1/M
        poisson_mat = np.matrix(np.zeros(C.shape))
        for i in range(C.shape[0]):
            for j in range(M):
                poisson_mat[i,j] = np.random.poisson(C[i,j])
        if norm: # seems to make no performance difference 1/25
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
        corr = np.nan_to_num(np.corrcoef(poisson_mat))
        return corr
    avg_result = (reduce(operator.add, (poisson_corr(mat, norm=norm) for i in
                                        range(repeat))) / repeat)
    return avg_result

