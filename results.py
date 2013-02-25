import numpy as np
import random
import utils as ut
from Struct import Struct
import ml
import features as fe
import cluster as cl
import ppi
import cyto
import compare as cp
import os
import pairdict as pd
import corum as co
import examples as ex

clf_factories = {'svm': (ml.svm, ml.linear), 'tree': (ml.tree, ml.tree_feats)}

def cvtest(name, base_sp, nsp, fs, base_featstruct, kfold=10, clf_type='tree',
        nfeats=40, norm=True, ppi_output=None, train_limit=None,
        save_data=True, balance_train=True, **kwargs):
    """
    
    """
    exs = ppi.feature_array(base_sp, fs, base_featstruct,
            nsp, **kwargs) if ppi_output is None else ppi_output
    arrfeats, ntest_pos = fe.arr_copy(exs.arrfeats), exs.ntest_pos
    assert len(arrfeats)>0, '0 examples not supported'
    if train_limit: print 'Sampling %s train/cv examples' % train_limit
    train_limit = train_limit or len(arrfeats)
    arrfeats = fe.keep_rows(arrfeats, random.sample(range(len(arrfeats)),
        int(train_limit))) # shuffle even if not sampling. don't random.shuffle
    ntest_pos = int(ntest_pos * train_limit / len(arrfeats))
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
            source_feats=exs.arrfeats.dtype.names)
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

def deprecated_test(name, base_sp, nsp, fs, ttbase, clf=None,
        clffact_feats=None, nfeats=40, norm=True, ppi_output=None,
        train_limit=None, save_data=True, **kwargs):
    """
    Making base train/test: give ttbase as None, pass do_filter=False, fs=[],
    extdata=[], nfeats=None
    """
    exs = ppi.learning_examples(base_sp, fs, ttbase,
            nsp, **kwargs) if ppi_output is None else ppi_output
    feats = feature_selection(exs.train, nfeats, clffact_feats)
    arr_train, arr_test = [fe.keep_cols(a,feats) for a in exs.train, exs.test]
    # Careful: clf is length 0 even when instantiated. Do not say 'if clf:'
    clf = clf if clf is not None else ml.svm()
    print "Classifier:", clf
    if train_limit: 
        print 'Sampling %s training examples' % train_limit
        arr_train = fe.keep_rows(arr_train,
                random.sample(range(len(arr_train)),train_limit))
    scaler = ml.fit_clf(arr_train, clf, norm=norm)
    if ml.exist_pos_neg(arr_train):
        ppis = ml.classify(clf, arr_test, scaler=scaler)
    else:
        ppis = []
    result = Struct(train=arr_train[['id1','id2','hit']], clf=clf,
            scaler=scaler, ppis=ppis, ntest_pos=exs.ntest_pos, name=name,
            species=base_sp, ppi_params=str(clf), feats=feats,
            source_feats=exs.train.dtype.names)
    if save_data:
        result.exs = exs
    return result

def feature_selection(arr, nfeats, clf):
    if nfeats>0:
        print 'Selecting top %s of %s features' % (nfeats,
            len(arr.dtype.names)-3)
        feats = ml.feature_selection(arr, clf)[0][:nfeats]
    else:
        feats = list(arr.dtype.names[3:])
    return feats

def predict(name, sp, arrsource, arrfeats, nsp, clf_scaler_feats=None,
        clf_base=None, clf_feats=None, norm=True, nfeats=40,
        balance_train=False):
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
        feats = feature_selection(arrfeats, nfeats, clf_feats)
        arrfeats = fe.keep_cols(arrfeats, feats)
        clf = clf_base if clf_base is not None else ml.svm()
        scaler = ml.fit_clf(arrfeats, clf, norm=norm)
    print "Classifier:", clf
    arrsource = fe.keep_cols(arrsource, feats)
    ppis = ml.classify(clf, arrsource, scaler=scaler)
    return Struct(ppis=ppis,name=name, species=sp, ppi_params=str(clf),
            feats=feats, nsp=nsp, arrfeats=arrfeats)

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

def predict_clust(name, sp, scored, exs, nsp, savef=None, pres=None, 
        pd_spcounts=None, **kwargs):
    savef = savef if savef else ut.bigd(name)+'.pyd'
    if pres is None:
        pres_fname = ut.pre_ext(savef, '_pres')
        assert not os.path.exists(pres_fname), "Destination filename exists"
        pres = predict(name, sp, scored, exs, nsp, **kwargs)
        pres.splits = exs.splits
        ut.savepy(pres, pres_fname) 
    clusts = cl.multi_clust(pres.ppis)
    clstruct = cp.result_stats(sp, exs.splits, clusts, nsp)
    ut.savepy(clstruct, ut.pre_ext(savef, '_clstruct'))
    pres.cxs, pres.cxppis, pres.ind = cp.select_best(clstruct, ['ppv','mmr'])
    ut.savepy([pres.cxs,pres.cxppis],
    ut.pre_ext(savef,'_cxs_cxppis_clust%s_%scxs' % (pres.ind, len(pres.cxs))))
    cyto_export(pres, pres.train, name_ext='_clust%s_%scxs' % (pres.ind,
        len(pres.cxs)), geneatts=ut.proj_path('gene_desc_'+sp),
        pd_spcounts=pd_spcounts)
    return pres, clstruct

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

def cyto_export(result, arrtrain, ppis=None, cxs='default', geneatts=
        'Hs_ensg_name_desc_uni_entrez.tab', species=None, pd_spcounts=None, 
        name_ext=''):
    fname = 'cy_'+result.name+name_ext+'.tab'
    species = species if species else result.species
    ppis = ppis if ppis else result.cxppis
    cxs = cxs if cxs!='default' else result.cxs
    cyto.cyto_prep(ppis, arrtrain, fname, geneatts, cxs, species=species,
            pd_spcounts=pd_spcounts)

def replace_eluts(arr, scores, name='eluts'):
    species = set(ut.config()['elut_species'].split('_'))
    keep_cols = [n for n in arr.dtype.names[3:] if n[:2] not in species]
    keep_cols = keep_cols + [name]
    newarr = ppi.base_array(arr[['id1','id2','hit']], keep_cols, len(arr))
    newarr[name] = scores
    return newarr
