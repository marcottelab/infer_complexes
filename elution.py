from __future__ import division
import itertools
import numpy as np
import operator
import os
import random
import ml
import orth
import pairdict
from pairdict import PairDict
import score
from Struct import Struct
import utils as ut

def load_elution(fname, getname=True):
    # expected file structure:
    # first col: gene id
    # second col: treat differently if 2nd col header is 'Total' or
    # 'Description'
    # remaining cols: elution profile data
    lines = [l for l in ut.load_tab_file(fname)]
    # final row: total count in msblender output; don't skip in cuihong's data
    skip_final_row = (lines[-1][0][0] == '#')
    rows = lines[1:-1] if skip_final_row else lines[1:]
    fractions = [f for f in lines[0][1:]]
    if fractions[0].lower() in ['total', 'totalcount', 'description']:
        start_data_col = 2
        fractions.remove(fractions[0])
    else:
        start_data_col = 1
    mat = np.matrix([row[start_data_col:] for row in rows],dtype='float64')
    prots = [row[0] for row in rows]
    elut = Struct(mat=mat, prots=prots, fractions=fractions, filename=fname,
                  filename_original=fname)
    if start_data_col == 2:
        col2name_vals = [row[1] for row in rows]
        elut.column2vals = col2name_vals
    if getname: elut.name = os.path.basename(fname).split('.')[0]
    return elut

def all_filtered_pairs(fnames, score_keys, cutoff=0.25, sp_base=None,
                       seqdb=None, verbose=True):
    allpairs = PairDict([])
    for skey,f in itertools.product(score_keys,fnames):
        if verbose: print skey, cutoff, f
        elut = load_elution(f)
        pair_inds = score.pairs_exceeding(elut, skey, thresh=cutoff)
        newpairs = PairDict(((elut.prots[i], elut.prots[j])
                             for (i,j) in pair_inds))
        newpairs = translate_pairs(newpairs, sp_base, file_sp(f), seqdb)
        allpairs = pairdict.pd_union_novals(allpairs, newpairs)
    return allpairs

def translate_pairs(pairs, sp_base, sp_target, seqdb):
    t2b = targ2base(sp_base, sp_target, seqdb)
    if t2b:
        pairs = PairDict(((base1,base2)
                          for t1,t2 in pairs.d
                          for base1,base2 in
                          itertools.product(t2b.get(t1,[]), t2b.get(t2,[]))))
    return pairs
    
def targ2base(sp_base, sp_target, seqdb):
    t2b = None
    if sp_base:
        if sp_base != sp_target:
            t2b = orth.odict(sp_target+'_'+seqdb, sp_base+'_'+seqdb)
    return t2b

def file_sp(filename):
    return ut.shortname(filename)[:2]

def all_prots(elut_fs, sp_base=None, seqdb=None):
    allprots = set([])
    for f in elut_fs:
        t2b = targ2base(sp_base, file_sp(f), seqdb)
        newprots = set(load_elution(f).prots) 
        if t2b:
            newprots = set((b for t in newprots for b in t2b.get(t,[])))
        allprots = set.union(allprots, newprots)
    return allprots


##################################################
# Special-purpose or single-use
##################################################

def _fraction_elutions(fractions):
    """
    Given a list of fraction names, group them after removing 'FractionXX' from
    the end.  Return a dict of { elutionname: listofindices }.
    Example of fraction name: Orbitrap_HeLaCE_IEF_pH3_to_10_Fraction10
    """
    elution_names = {}
    for i,fname in enumerate(fractions):
        ename = fname[:fname.find('_Fraction')]
        elution_names.setdefault(ename,[]).append(i)
    return elution_names
    
def split_muliple_elutions(big_elut):
    """
    Split an elution into multiple based on use of _fraction_elutions.
    """
    elution_columns = _fraction_elutions(big_elut.fractions)
    eluts = {}
    for elution_name in elution_columns:
        new_elut = Struct()
        new_elut.__dict__ = big_elut.__dict__.copy()
        new_elut.mat = big_elut.mat[:,elution_columns[elution_name]]
        new_elut.fractions = list(np.array(big_elut.fractions)[elution_columns[elution_name]])
        new_elut.filename = big_elut.filename + '__' + elution_name
        eluts[elution_name] = new_elut
    return eluts

def write_elution(elut, fname, forR=False):
    """
    Write out an elution in the spcount format
    $ProtID\tTotalCount\tCol1....
    """
    # First eliminate empty protein rows
    nonzeros = np.sum(np.array(elut.mat),axis=1)>0
    arr = np.array(elut.mat[nonzeros,:])
    prots = list(np.array(elut.prots)[nonzeros])
    if not forR:
        header = "#ProtID TotalCount".split() + elut.fractions
        data = [[prots[i], np.sum(mat[i,:])] + arr[i,:].tolist() for i in
                range(len(prots))]
    else: #R: no column header for first column, and transpose
        header = prots
        data = [[elut.fractions[i]] + arr[:,i].tolist() for i in
                range(len(elut.fractions))]
    ut.write_tab_file([header] + data, fname)

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

def downsample_elution(elution, downsample, seed=0):
    """
    Return a new elution with every downsample-th fraction.
    """
    down_elut = Struct()
    down_elut.__dict__ = elution.__dict__.copy()
    down_elut.mat = elution.mat[:,seed::2]
    down_elut.fractions = elution.fractions[::2]
    down_elut.name = elution.name + '_down%i' % downsample
    return(down_elut)
