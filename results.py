from __future__ import division
import itertools as it
import numpy as np
import os
import random

import cluster as cl
import cyto
import compare as cp
import corum as co
import cv
import evidence as ev
import examples as ex
import features as fe
import ml
import pairdict as pd
import ppi
from Struct import Struct
import utils as ut

clf_factories = {'svm': (ml.svm, ml.linear), 'tree': (ml.tree, ml.tree_feats)}

def cvtest(name, base_sp, nsp, fs, base_featstruct, kfold=2, clf_type='svm',
        nfeats=100, norm=True, ppi_output=None, train_limit=None,
        save_data=True, balance_train=False, keep_cols=None, clf_factory=None,
        clffact_feats=None, **kwargs):
    """
    """
    assert kfold>1, "CV K-fold 1 not possible"
    exs = ppi.feature_array(base_sp, fs, base_featstruct,
            nsp, **kwargs) if ppi_output is None else ppi_output
    arrfeats, ntest_pos = fe.arr_copy(exs.arrfeats), exs.ntest_pos
    assert len(arrfeats)>0, '0 examples not supported'
    if train_limit: print 'Sampling %s train/cv examples' % train_limit
    train_limit = train_limit or len(arrfeats)
    arrfeats = arrfeats if keep_cols is None else fe.keep_cols(arrfeats, keep_cols)
    arrfeats = fe.keep_rows(arrfeats, random.sample(range(len(arrfeats)),
        int(train_limit))) # shuffle even if not sampling. don't random.shuffle
    ntest_pos = int(ntest_pos * train_limit / len(arrfeats))
    if clf_type in clf_factories and clf_factory is None:
        clf_factory, clffact_feats = clf_factories[clf_type]
    ppis = []
    for k in range(kfold):
        print 'Fold %s:' % k
        ppis_fold,clf,scaler,feats = fold_test(arrfeats, kfold, k, clf_factory,
                clffact_feats, nfeats, norm, balance_train)
        ppis += ppis_fold
    random.shuffle(ppis)
    ppis.sort(key=lambda x: x[2], reverse=True)
    result = Struct(traincv=arrfeats[['id1','id2','hit']], clf=clf,
            scaler=scaler, ppis=ppis, ntest_pos=ntest_pos, name=name,
            species=base_sp, ppi_params=str(clf), feats=feats,
            source_feats=exs.arrfeats.dtype.names, balance_train=balance_train)
    if save_data:
        result.exs = exs
    return result


def fold_test(arrfeats, kfold, k, clf_factory, clffact_feats, nfeats, norm,
        balance_train):
    arrtrain, arrtest = fe.arr_kfold(arrfeats, kfold, k)
    if balance_train:
        arrtrain = fe.balance_train(arrtrain)
    if nfeats:
        clf_feats = clffact_feats()
        feats = feature_selection(arrtrain, nfeats, clf_feats)
        arrtrain,arrtest = [fe.keep_cols(a,feats) for a in arrtrain,arrtest]
    else:
        feats = None
    clf = clf_factory()
    if k==0: print "Classifier:", clf
    scaler = ml.fit_clf(arrtrain, clf, norm=norm)
    if ml.exist_pos_neg(arrtrain):
        ppis = ml.classify(clf, arrtest, scaler=scaler, do_sort=False)
    else:
        ppis = []
    return ppis,clf,scaler,feats

def feature_selection(arr, nfeats, clf):
    if clf is not None and nfeats>0:
        print 'Selecting top %s of %s features' % (nfeats,
            len(arr.dtype.names)-3)
        feats = ml.feature_selection(arr, clf)[0][:nfeats]
    else:
        feats = list(arr.dtype.names[3:])
    return feats

def predict(name, sp, arrsource, arrfeats, nsp, clf_scaler_feats=None,
        clf_factory=None, clffact_feats=None, clf_type='svm', norm=True, nfeats=100, balance_train=False):
    """
    - arrfeats: labeled training examples array, from
      ppi.feature_array.arrfeats, also stored in res.cvtest result as
      result.exs.arrfeats.
    - arrsource: array of data to classify, matching the training array
    """
    if clf_scaler_feats:
        clf, scaler, feats = clf_scaler_feats
    else:
        if balance_train:
            arrfeats = fe.balance_train(arrfeats)
        if clf_type in clf_factories and clf_factory is None:
            clf_factory, clffact_feats = clf_factories[clf_type]
        feats = feature_selection(arrfeats, nfeats, clffact_feats() if clffact_feats
                else None)
        arrfeats = fe.keep_cols(arrfeats, feats)
        clf = clf_factory()
        scaler = ml.fit_clf(arrfeats, clf, norm=norm)
    print "Classifier:", clf
    arrsource = fe.keep_cols(arrsource, feats)
    ppis = ml.classify(clf, arrsource, scaler=scaler)
    pres = Struct(ppis=ppis,name=name, species=sp, ppi_params=str(clf),
            feats=feats, nsp=nsp, arrfeats=arrfeats,
            balance_train=balance_train)
    return pres

def combine_cvres_ppis(cresa, cresb):
    """
    ntest_pos is fixed and updated.  splits are not carried over.
    """
    newres = combine_pres_ppis(cresa, cresb)
    lpa,lpb = [res.exs.splits[0][0] for res in cresa,cresb]
    merged_lpairset = ex.merge_lps([lpa,lpb])
    newres.ntest_pos = len(merged_lpairset.pairs)
    return newres

def combine_pres_ppis(resa, resb):
    res = ut.struct_copy(resa)
    res.name = 'combined: %s, %s' % (resa.name, resb.name)
    res.ppis = pd.pd_lol(pd.pd_combine_ppis(pd.PairDict(resa.ppis),
            pd.PairDict(resb.ppis), combine_or))
    res.ppis.sort(key=lambda x: x[2], reverse=True)
    return res

def combine_ppis(a,b):
    ppis = pd.pd_lol(pd.pd_combine_ppis(pd.PairDict(a), pd.PairDict(b),
        combine_or))
    return ppis

def combine_or(a,b):
    return 1-((1-a)*(1-b))

def combine_ppis_matched(ppisa, ppisb):
    return [(pa[0],pa[1],combine_or(pa[2],pb[2]),pa[3])
            for pa,pb in ut.zip_exact(ppisa, ppisb)]

def predict_clust(name, sp, nsp, obs=None, exs=None, savef=None, pres=None,
        pd_spcounts=None, cl_kwargs={}, clusts=None, runid=None,
        count_ext=False, cutoff=0.5, n_cvs=7, accept_clust=False,
        obs_fnames=None, base_splits=None, obs_kwargs={}, kfold=3,
        gold_nspecies=2, do_cluster=True, do_2stage_cluster=True,
        cxs_cxppis=None, do_rescue=True, n_rescue=20000, rescue_fracs=20,
        rescue_score=0.9, clstruct=None, **predict_kwargs):
    """
    - obs/test_kwargs: note obs_kwargs is combined with predict_kwargs to enforce
      consistency.
    - pd_spcounts: supply from ppi.predict_all if nsp > 1.
    - base_splits: supply exs.splits to generate examples from existing
      division of complexes.
    - cxs_cxppis: provide if you want to export, or do the ppi rescue
      clustering--also must set accept_clust=True, do_rescue=True
    """
    savef = savef if savef else ut.bigd(name)+'.pyd'
    print "Will save output to", savef
    runid = runid or random.randrange(0,1000)
    if clusts is None: 
        if pres is None:
            if obs is None:
                obs, pd_spcounts = ppi.predict_all(sp, obs_fnames,
                        save_fname=savef.replace('.pyd',''), nsp=nsp,
                        **obs_kwargs)
            if exs is None:
                cvtest_kwargs = ut.dict_quick_merge(obs_kwargs, predict_kwargs)
                n_cvs = 1 if base_splits is not None else n_cvs
                cvs, cvstd = cvstd_via_median(name, sp, nsp, obs_fnames, kfold,
                        base_splits, n_cvs, **cvtest_kwargs)
                if n_cvs > 1:
                    ut.savepy(cvs, ut.pre_ext(savef, '_cvs_%s' % n_cvs))
                ut.savepy(cvstd, ut.pre_ext(savef, '_cvstd'))
                exs=cvstd.exs
            pres = predict(name, sp, obs, exs.arrfeats, nsp, **predict_kwargs)
            pres.exs = exs
            ut.savepy(pres, ut.pre_ext(savef, '_pres'), check_exists=True) 
        else:
            pres=ut.struct_copy(pres)
            if do_rescue:
                assert obs is not None, "Must supply obs for rescue step"
    merged_splits = pres.exs.splits[1] # splits is (lp_splits, clean_splits)
    if do_cluster:
        if cxs_cxppis is None and clstruct is None:
            if clusts is None and cxs_cxppis is None:
                #if calc_fracs:
                    #cl_kwargs['fracs'] = [cp.find_inflection(pres.ppis, merged_splits,
                        #pres.species, gold_nspecies)]
                clusts, runid = multi_clust(pres.ppis, savef=savef, runid=runid,
                        pres=pres, gold_splits=merged_splits,
                        gold_nspecies=gold_nspecies, **cl_kwargs)
                ut.savepy(clusts, ut.pre_ext(savef, '_clusts_id%s' % runid))
            if do_2stage_cluster:
                clusts2 = multi_stage2_clust(clusts, pres.ppis, runid=runid,
                        **cl_kwargs)
                clstruct = cp.result_stats(sp, merged_splits, clusts2,
                        gold_nspecies) 
                ut.savepy(clstruct, ut.pre_ext(savef, '_clstruct2_id%s' % runid))
            else:
                clstruct = cp.result_stats(sp, merged_splits, clusts, nsp) 
                ut.savepy(clstruct, ut.pre_ext(savef, '_clstruct_id%s' % runid))
        if accept_clust:
            if cxs_cxppis is None:
                pres.cxs, pres.cxppis, pres.ind = cp.select_best(clstruct)
                ut.savepy([pres.cxs,pres.cxppis],
                        ut.pre_ext(savef,'_cxs_cxppis_id%s_ind%s_%scxs'
                            % (runid, pres.ind, len(pres.cxs))))
            else:
                pres.cxs, pres.cxppis = cxs_cxppis
                pres.ind = 0
            if do_rescue:
                # note cl_kwargs aren't passed--would be messy
                pres.cxs, pres.cxppis, pres.ppis_rescue = rescue_ppis(pres,
                        obs, n_rescue, cutoff_fracs=rescue_fracs,
                        cutoff_score=rescue_score)
            cyto_export(pres, merged_splits, name_ext='_clust%s_%scxs' % (pres.ind,
                len(pres.cxs)), pd_spcounts=pd_spcounts, arrdata=obs,
                cutoff=cutoff, count_ext=False, arrdata_ppis=None)
            return pres
        else:
            return pres, clstruct
    else:
        return pres

def rescue_ppis(pres, obs, n_rescue, cutoff_fracs=20, cutoff_score=0.9,
        exclude_ppis=None, frac_retain=0.1, cltype='mcl', merge_cutoff=0.55,
        I=3, negmult=1):
    ppis_counted = ev.ppis_fracspassing_counts(pres.ppis[:n_rescue], obs,
            exclude_ppis=(exclude_ppis or pres.cxppis))
    print "%s PPIs considered for possible rescue" % len(ppis_counted)
    ppis_rescue = [p for p in ppis_counted 
            if p[-1] > cutoff_fracs or p[2] > cutoff_score]
    print "%s ppis exceeding score cutoff, %s exceeding fracs cutoff" % (len([p
        for p in ppis_counted if p[-1] > cutoff_fracs]), len([p for p in
            ppis_counted if p[2] > cutoff_score]))
    print "%s PPIs to be rescued" % len(ppis_rescue)
    ppis_rescue_as_cxs = [set((p[0],p[1])) for p in ppis_rescue]
    additional_cxs = pres.cxs + ppis_rescue_as_cxs
    cxstruct_rescue = cl.filter_clust(ppis_rescue, ut.list_frac(pres.ppis,
        frac_retain), cltype='mcl', merge_cutoff=merge_cutoff, I=I,
        negmult=negmult, add_to_cxs=additional_cxs)
    print "Total complexes >= size 3: was %s, now %s" % (len([c for c in
        pres.cxs if len(c)>=3]), len([c for c in cxstruct_rescue.cxs if len(c)>=3]))
    return cxstruct_rescue.cxs, cxstruct_rescue.cxppis, ppis_rescue

def cvstd_via_median(name, sp, nsp, obs_fnames, kfold, base_splits, n_cvs,
        savef=None, overwrite_balanced=True, **kwargs):
    """
    overwrite_balanced: If true, overwrite balance_train parameter to be true
    to speed up this selection process, as it doesn't matter.
    """
    if overwrite_balanced:
        kwargs['balance_train'] = True
    cvs = [cvtest(name, sp, nsp, obs_fnames, None, kfold=kfold,
        both_cx_splits=base_splits, **kwargs) for i in
        range(n_cvs)]
    if base_splits is not None:
        cvstd = cvs[0] # only one, built from provided base_splits
    else:
        cvstd = median_cvtest(cvs)
    return cvs, cvstd

def median_cvtest(cv_results):
    auprs = [(i, cv.aupr(c.ppis, c.ntest_pos)) for i,c in
            enumerate(cv_results)]
    auprs.sort(key=lambda x: x[1])
    median_index = auprs[int(len(auprs)/2)][0]
    return cv_results[median_index]

def multi_stage2_clust(clusts, ppis_all, runid=None, frac_retain=.1,
        I_params=[2,3,4], **kwargs):
    clusts2 = []
    runid = runid or random.randrange(1,1000)
    for cxstruct1 in clusts:
        for I in I_params:
            ppis_retain = ut.list_frac(ppis_all, frac_retain)
            cxstruct2 = cl.stage2_clust(cxstruct1.cxppis, ppis_retain,
                    cxstruct1.cxs, runid=runid, I=I, **kwargs)
            cxstruct2.params = cxstruct1.params + ",I_mcl=%s" % I
            clusts2.append(cxstruct2)
    return clusts2

def multi_clust(tested, score_cutoffs=None, length_cutoffs=None,
        fracs=[.012,.014], frac_retain=.1, ds=[.1,.25,.3,.35], ms=[.1,.15,.2],
        penalties=[.1,1], overlaps=[.55], haircuts=[0,.2], max_pval=1,
        savef=None, runid=None, show_stats=True, pres=None, gold_nspecies=1,
        gold_splits=None, gold_minlen=3, mdprod_min=.01, **kwargs):
    runid = runid or random.randrange(1,1000)
    fracs = (fracs if fracs is not None 
        else [cl.n_thresh(tested, s)/len(tested) for s in score_cutoffs] if score_cutoffs is not None
        else [le/len(tested) for le in length_cutoffs])
    print "random id:", runid
    clusts = []
    params = [fracs, ds, ms, penalties, overlaps, haircuts]
    products = it.product(*params)
    for (f,d,m,p,o,h) in products:
        if d*m >= mdprod_min:
            cxstruct = cl.filter_clust(ut.list_frac(tested, f),
                    ut.list_frac(tested, frac_retain), merge_cutoff=o, negmult=m, min_density=d,
                    runid=runid, penalty=p, max_pval=max_pval, max_overlap=o,
                    haircut=h, **kwargs)
            cxstruct.params = ('density=%s,frac=%s,f_retain=%s,negmult=%s,penalty=%s,max_overlap=%s,haircut=%s' % (d,f,frac_retain,m,p,o,h))
            clusts.append(cxstruct)
            if show_stats and len(cxstruct.cxs)>0:
                if pres is not None and gold_splits is not None:
                    out = cp.select_best(cp.result_stats(pres.species, gold_splits,
                        clusts[-1:], gold_nspecies, min_gold_size=gold_minlen))
                else:
                    print "Can't show stats: pres and gold_splits required."
            if savef and (len(clusts) % 10 == 1):
                ut.savepy(clusts, ut.pre_ext(savef, "clusts_temp_%s_%s" % (ut.date(),
                    runid)))
    return clusts, runid

def combine_train(atrain, atest):
    #ltrain, ltest = [[r for r in arr if r[2]==1] for arr in atrain,atest]
    inds_test_trues = np.where(atest['hit']==1)[0]
    inds_test_falses = np.where(atest['hit']==0)[0]
    if len(inds_test_trues) < len(inds_test_falses): 
        inds_test_falses = random.sample(inds_test_falses, len(inds_test_trues)) 
    inds_all = np.concatenate((inds_test_trues, inds_test_falses))
    inds_all.sort()
    comb = np.concatenate((atrain, atest[inds_all]))
    return comb

def multi_base_tests(sp, nsp, fs, ntests=5, **kwargs):
    reslist = []
    for i in range(ntests): 
        exs = ppi.learning_examples(sp,[],None,nsp,do_filter=False,extdata=[])
        restest = test('%s %ssp, %s of %s' %(sp,nsp,i,ntests),sp,nsp,fs,exs,
                **kwargs)
        reslist.append((exs,restest))
    return reslist

def plot_reslist(l_bases_tests):
    import plotting as pl
    for base,rtest in l_bases_tests:
        pl.plot_result(rtest)

def cluster(result, fracppis, **kwargs):
    """
    Works on output of either test or predict.
    """
    nppis = int(fracppis*len(result.ppis))
    result.cxs, result.cxppis = cl.filter_c1(result.ppis[:nppis], **kwargs)
    result.cx_params = str(kwargs)
    return result

def cyto_export(result, cxs_splits, ppis=None, cxs='default', species=None,
        pd_spcounts=None, name_ext='', arrdata=None, do_splist=True,
        cutoff=0.5, count_ext=False, **kwargs):
    fname = 'cy_'+result.name+name_ext+'.tab'
    species = species if species else result.species
    ppis = ppis if ppis else result.cxppis
    cxs = cxs if cxs!='default' else result.cxs
    cyto.cyto_prep(ppis, cxs_splits, fname, cxs, species=species,
            pd_spcounts=pd_spcounts, do_splist=do_splist, arrdata=arrdata,
            **kwargs)

def replace_eluts(arr, scores, name='eluts'):
    species = set(ut.config()['elut_species'].split('_'))
    keep_cols = [n for n in arr.dtype.names[3:] if n[:2] not in species]
    keep_cols = keep_cols + [name]
    newarr = ppi.base_array(arr[['id1','id2','hit']], keep_cols, len(arr))
    newarr[name] = scores
    return newarr

