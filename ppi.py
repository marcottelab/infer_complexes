from __future__ import division
import examples as ex
import corum as co
import seqs
import orth
import utils as ut
import elution as el
from Struct import Struct
import numpy as np
import fnet
import score
from numpy import zeros,ndarray
from pairdict import PairDict
import features as fe
import os
extfs = ['ext_Dm_guru','ext_Hs_malo'] 

def feature_array(species, elut_fs, base_exs, nsp,
        scores=['poisson','wcc','apex','pq_euc'], extdata=['net_Hs19']+extfs,
        split_props=[0,.5,1], ind_cycle=None, cutoff=0.5, both_cx_splits=None,
        gold_consv_sp='Dm', do_filter=True, gidscheme='', go_location=None,
        filter_incl_ext=False, postfilt_sps=None, allow_singles='auto',
        do_filter_corum=True, do_post_filter_override=False, **score_kwargs):
    """
    - Species: 'Hs' or 'Ce'... The base of the predictions in terms of
      identifiers used and orthology pairings.
    - base_exs: uses test, train, splits from output of previous run.
    - extdata: list of items that are either a string data_key matching to
      those specified in config.py, or a tuple ('name', data).
      The data should simply be a sequence of list/tuple of (id1, id2, score).
      other_evidence just appends to this.
    - do_filter: only set to false if you want a full unfiltered test/train arr
    - require_base: if True, only keeps interactions with score exceeding
      cutoff in the base species.  if False, in any species.
    - ind_cycle: leave as none for stochastic rotation when assigning positives
      to splits; set as eg [0,1] for cycling repeatedly through 2 splits.
    - Making base train/test: give ttbase as None, pass do_filter=False, fs=[],
      extdata=[].
    - Postfilt_sps: after everything, retain only columns from the provided
      list of eg ['Hs','HS']. Removes empty rows after, but no other filtering after.
      """
    gold_consv_sp = gold_consv_sp if nsp>1 else '' #ignore for 1-sp case
    if nsp > 1 and do_filter and not filter_incl_ext: 
        check_nspecies(elut_fs, nsp)
    if allow_singles=='auto':
        allow_singles = (nsp > 1)
    if base_exs:
        pdtrain = (pd_from_arr(base_exs.arrfeats) if
                hasattr(base_exs,'arrfeats') else base_exs.pdtrain)
        splits = base_exs.splits
        ntest_pos = base_exs.ntest_pos
    else:
        print "\n\nGenerating training examples from %s splits." % ("NEW" if
        both_cx_splits is None else "PROVIDED")
        pdtrain,splits = base_splits(species, elut_fs, split_props,
                ind_cycle, both_cx_splits, gold_consv_sp, gidscheme,
                go_location, do_filter_corum)
        ntest_pos = len([v for k,v in pdtrain.d.items() if v[0]==1])
    print 'total train/cv positives:', ntest_pos
    arrfeats = new_score_array(pdtrain, scores, elut_fs, extdata) 
    arrfeats = score_and_filter(arrfeats, scores, elut_fs, cutoff, species,
            extdata, gidscheme, nsp=nsp, do_filter=do_filter,
            filter_incl_ext=filter_incl_ext, postfilt_sps=postfilt_sps,
            allow_singles=allow_singles, **score_kwargs)
    if nsp > 1 and (do_filter or do_post_filter_override): # why do this after the previous step? 3/8/13
        arrfeats = fe.filter_nsp_nocounts(arrfeats, nsp=nsp, cutoff=cutoff,
                count_ext=filter_incl_ext) 
    print 'done.', stats(arrfeats)
    return Struct(arrfeats=arrfeats, ntest_pos=ntest_pos, splits=splits)


def check_nspecies(fnames, nsp):
    nspec_fracs = len(set([ut.shortname(f)[:2] for f in fnames]))
    assert nspec_fracs >= nsp, ("** Can't use %s species; only %s species" %
            (nsp,nspec_fracs))

def pd_from_arr(arr):
    return PairDict([[p[0],p[1],p[2]] for p in arr])

def base_splits(species, elut_fs, splits, ind_cycle, both_cx_splits, consv_sp,
        gidscheme, go_location, do_filter_corum):
    ppi_cxs,clean_cxs,all_cxs = load_training_complexes(species,
            gidscheme, consv_sp=consv_sp, go_location=go_location,
            do_filter_methods=do_filter_corum)
    return ex.base_examples_single(ppi_cxs, clean_cxs, all_cxs,
            splits, both_cx_splits=both_cx_splits, ind_cycle=ind_cycle)

def predict_all(species, elut_fs, scores=['poisson','wcc','apex','pq_euc'],
        extdata=['net_Hs19']+extfs, nsp=1, cutoff=0.5, base_arr=None,
        do_filter_sp=True, save_fname=None, gidscheme='', postfilt_sps=None,
        filter_incl_ext=False, verbose=False, allow_singles='auto',
        **score_kwargs):
    """
    Same more or less as learning_examples above, but produces all predictions
    in the elution files.
    pd_spcounts will be None if nsp==1.
    """
    if allow_singles=='auto':
        allow_singles = (nsp > 1)
    if not base_arr:
        pd_all = el.all_filtered_pairs(elut_fs, scores, cutoff==cutoff,
                sp_base=species, allow_singles=allow_singles)
        print len(pd_all.d), 'total interactions passing cutoff'
        for k in pd_all.d: pd_all.d[k] = -1 #interaction marked as unknown
        arr = new_score_array(pd_all, scores, elut_fs, extdata)
        del pd_all #lots of memory
    else:
        arr = base_arr
    scored_arr = score_and_filter(arr, scores, elut_fs, cutoff, species,
                extdata, gidscheme, nsp=nsp, filter_incl_ext=filter_incl_ext,
                postfilt_sps=postfilt_sps, allow_singles=allow_singles,
                **score_kwargs)
    if save_fname: np.save(save_fname+'1sp.npy', scored_arr)
    spcounts = None
    if nsp > 1 and do_filter_sp:
        scored_arr, spcounts = fe.filter_nsp(scored_arr, nsp=nsp,cutoff=cutoff,
                count_ext=filter_incl_ext)
        if save_fname: 
            np.save(save_fname+str(nsp)+'sp.npy', scored_arr)
            ut.savepy(spcounts, save_fname+'_pdspcounts%s.pyd'%cutoff)
    if verbose:
        print scored_arr.dtype.names
    return scored_arr, spcounts

def score_and_filter(arr, scores, elut_fs, cutoff, species, extdata, gidscheme,
        do_filter=True, require_base=None, single_base_orth=False,
        filter_multi_orths=0.25, nsp=1, filter_incl_ext=False,
        postfilt_sps=None, allow_singles=True):
    require_base = require_base if require_base is not None else (nsp==1)
    filter_multi_orths = (filter_multi_orths if (nsp>1 and not require_base)
            else False) #ignore if 1sp
    print '\nScoring %s elutions with %s base, scores: %s.' % (len(elut_fs),
            species, ','.join(scores))
    # FIX: currently not passing gidscheme, as I can't handle that currently
    score.score_array_multi(arr, species, elut_fs, scores, cutoff,
            remove_multi_base=single_base_orth, gidscheme='')
    if do_filter:
        if cutoff != -1 and do_filter:
            print 'Filtering, require_base =', require_base
            require_species = set([species]) if require_base else None
            arr = fe.filter_require_sp(arr, require_species, cutoff=cutoff,
                    count_ext=filter_incl_ext)
        if filter_multi_orths: 
            arr = fe.filter_multi_orths(arr, species, filter_multi_orths)
    if extdata:
        print '\nScoring with external data:', ','.join(extdata)
        sp_idscheme = (species + '_' + gidscheme) if gidscheme else species
        for d in extdata: fnet.score_arr_ext(arr, sp_idscheme, d)
        if postfilt_sps:
            arr = fe.keep_cols(arr, fe.species_cols(arr, set(postfilt_sps),
                count_ext=True))
    if cutoff == -1 and do_filter: # filter this if we haven't before.
        arr = remove_empties(arr)
    return arr 


def filter_arr(arr, columns, cutoff):
    """
    Uses max: return a new array of all rows from input arr with at least one
    value in columns exceeding cutoff.
    """
    newarr = arr[[row for row in range(len(arr)) if max([val for col,val in
                enumerate(arr[row]) if col in columns]) > cutoff]]
    del arr
    return newarr

def new_score_array(pd, scores, fs, extdata):
    """
    Pre-allocate an array sized based on the elution files and extdata.
    arr[0] returns the first row of data.
    row[name_score(f,s)] returns the proper score
    row['id1','id2','hit'] index to the first three columns
    """
    data_names = [score.name_score(f,s) for f in fs for s in scores]
    for key in extdata:
        ext_names = fnet.fnet_names(ut.config()[key])
        if ext_names:
            data_names += ext_names
        else:
            data_names.append(key)
    arr = base_array(((p1,p2,lhit) for (p1,p2),lhit in pd.d.items()),
            data_names, len(pd.d))
    return arr

def base_array(triple, data_names, data_len):
    """
    triple should be a sequence of (id1, id2, hit)
    """
    arr = zeros(shape = data_len, dtype = ','.join(['a20']*2 + ['i1'] +
        ['f2']*len(data_names)))
    names = ['id1','id2','hit'] + data_names
    arr.dtype.names = tuple(names)
    for i,(p1,p2,hit) in enumerate(triple):
        row = arr[i]
        row['id1'],row['id2'],row['hit'] = p1,p2,hit
    return arr

def check_overlaps(pdtrain, pdtest):
    overlaps = len([p for p,v in pdtrain.d.items() if pdtest.contains(p)])
    if overlaps:
        print '**Common pairs between train/test:', overlaps
    else:
        print "No common pairs between train/test."

def stats(arr):
    nums =  [sum(arr['hit']==tf) for tf in [1,0]]
    return 'arrfeats: %sP/%sN' % tuple(nums)

def load_training_complexes(species, gidscheme, consv_sp,
        do_filter_methods=True, go_location=None):
    """
    Uses:
    - ppi: for ppi training purposes. all unmerged corum complexes, but with
      duplicates removed, and when there's a name conflict, only the complex
      with the most members is retained.
    - merged: for clustering selection. corum complexes are merged based on
      redundancy criteria.
    - all: for visualizing what's known at the end. no merging, no name bash
      allowed--all complexes are returned that pass the filters. may contain
      complete duplicates.
    """
    ppi_cxs = co.load_ppi_cxs(go_location=go_location,
            do_filter_methods=do_filter_methods)
    merged_cxs = co.load_merged_cxs(go_location=go_location,
            do_filter_methods=do_filter_methods)
    all_cxs = co.load_ppi_cxs(minlen=2,maxlen=1000, go_location=go_location ,
            do_filter_methods=do_filter_methods)
    ppi_cxs = co.keep_longest(ppi_cxs)
    ppi_cxs,merged_cxs,all_cxs = [convert_complexes(cxs, species, gidscheme, consv_sp)
            for cxs in [ppi_cxs,merged_cxs,all_cxs]]
    return ppi_cxs, merged_cxs, all_cxs

def convert_complexes(cxs, species, gidscheme, consv_sp='', do_dedupe=True):
    """
    Convert from uniprot to human default (unless already human uniprot), then
    to other base species if specified.
    If consv_sp (conserved species) is provided, restrict gold standard to only
    include interactions for proteins with orthologs in that species.
    """
    #co.print_rib_count(cxs, '0')
    if gidscheme != 'uni':
        dconv = ut.load_dict_sets(ut.proj_path('corum_convert'))
        ps = seqs.load_prots_from_fasta(ut.proj_path('fastadir', 'Hs.fasta'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    #co.print_rib_count(cxs, '1')
    if species!='Hs':
        dconv = orth.odict('Hs', species)
        ps = seqs.load_prots_from_fasta(ut.proj_path('fastadir', 
                                                    species+'.fasta'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    if consv_sp:
        print "Requiring complexes conserved in %s" % consv_sp
        conserved = orth.odict(species, consv_sp)
        cxs = [(c,set([g for g in gs if g in conserved])) for c,gs in cxs]
    cxs = [(name,ps) for name,ps in cxs if len(ps)>1]
    co.print_rib_count(cxs, '2')
    return cxs

def preds_thresh(tested, thresh):
    limit = None
    for i,t in enumerate(tested):
        if t[2]<thresh:
            print 'limit:', i, t
            limit = i
            break
    return limit

def remove_empties(arr):
    if len(arr)==2: # assume this is test-train
        return [remove_empties(a) for a in arr]
    newarr = arr[list(arr.dtype.names[3:])]
    return arr[np.where([max(newarr[i]) for i in range(len(newarr))])]
