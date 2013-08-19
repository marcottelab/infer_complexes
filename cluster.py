from __future__ import division
import numpy as np
import os
import operator
import random
import subprocess

import corum as co
import compare as cp
import cyto
import pairdict as pd
from Struct import Struct
import utils as ut

mcl_command = '%(mclpath)s %(fin)s -I %(I)s -o %(fout)s --abc'
# Note -F csv means for clusterOne we get detailed output
c1_command = 'java -jar %(c1path)s -s %(min_size)s -d %(min_density)s --haircut %(haircut)s --penalty %(penalty)s --seed-method %(seed_method)s -F csv --max-overlap %(max_overlap)s %(fluff)s %(fin)s > %(fout)s ' 

def filter_clust(ppis_cluster, ppis_retain, negmult=50, cltype='c1',
        merge_cutoff=0.55, min_cx_length=2, is_recluster=False, **kwargs):
    print "Negatives rescale factor:", negmult
    cxstruct = cluster(ppis_cluster, negmult, cltype, **kwargs)
    if is_recluster:
        cxstruct.numbered_cxs = cxstruct.cxs
        cxstruct.cxs = clean_numbered_cxs(cxstruct.numbered_cxs, min_cx_length)
    cxstruct.cxs = remove_dupes(cxstruct.cxs, merge_cutoff, npasses=1)
    # be sure to get rid of complete duplicates
    cxstruct.cxs = remove_dupes(cxstruct.cxs, 1, npasses=3) 
    cxstruct.cxs = remove_subcomplexes(cxstruct.cxs, max_size=3)
    cxstruct.cxppis = _filter_ints(ppis_retain, cxstruct.cxs)
    print "%s Deduped complexes; %s cxppis" % (len(cxstruct.cxs), len(cxstruct.cxppis))
    return cxstruct

def stage2_clust(ppis_cluster, ppis_retain, cxs, **kwargs):
    cxppis = _filter_ints(ppis_cluster, cxs)
    cxlabeled_ppis = [p[:4] for p in cyto.ppis_as_cxs(cxppis, cxs)]
    # Note the labeling/unlabeling of ppis_retain is handled in filter_clust
    # since is_recluster is True.
    cxstruct = filter_clust(cxlabeled_ppis, ppis_retain, cltype='mcl',
            is_recluster=True, merge_cutoff=.55, **kwargs)
    return cxstruct

def clean_numbered_cxs(numbered_cxs, min_cx_length):
    return [set([g.split('_')[1] for g in c]) 
            for c in numbered_cxs if len(c) >= min_cx_length]

def clean_numbered_ppis(numbered_ppis):
    return [set([g.split('_')[1] for g in c]) 
            for c in numbered_cxs if len(c) >= min_cx_length]
    
def _filter_ints(inlist, cxs):
    d_cxs = dict([(i,set(c)) for i,c in enumerate(cxs)]) 
    pairs = co.pairs_from_complexes(d_cxs)
    pdp = pd.PairDict(pairs)
    return [tup for tup in inlist if pdp.contains((tup[0],tup[1]))]

def combine_clstructs(a,b):
    newclusts = a.cxstructs + b.cxstructs
    newstats = np.concatenate((a.stats, b.stats))
    return Struct(cxstructs = newclusts, stats = newstats)

clust_defaults = {
    'I': 1.8,
    'min_size': 2,
    'min_density': 0.2,
    'haircut': 0.0,
    'penalty': 2.0,
    'seed_method': 'cliques',
    'max_overlap': 0.7,
    'fluff': '--no-fluff',
    'c1path': os.path.expanduser('~/Dropbox/complex/tools')+'/cluster_one.jar',
    'mclpath': os.path.expanduser('~/local/bin')+'/mcl',
    'fin':
        os.path.expanduser('/data/complex/predictions')+
        '/temp_interactions.tab',
    'fout': os.path.expanduser('/data/complex/predictions') +
        '/temp_complexes.csv'
    }

def thresh(tested, score):
    return tested[:n_thresh(tested, score)]

def n_thresh(tested, score):
    for i,t in enumerate(tested):
        if t[2] < score:
            return i
    
def cluster(tested, negmult, cltype, max_pval=0.2, **kwargs):
    kwargs = ut.dict_set_defaults(kwargs, clust_defaults)
    command = c1_command if cltype=='c1' else mcl_command
    if 'runid' in kwargs: # keep temp files separate
        runid = str(kwargs['runid']) 
        kwargs['fin'] = ut.pre_ext(kwargs['fin'], runid)
        kwargs['fout'] = ut.pre_ext(kwargs['fout'], runid)
    export_cxs(tested, kwargs['fin'],negmult)
    command = command % kwargs
    print command
    shell_call(command)
    if cltype=='c1':
        cxs, pvals, cx_details = read_clust_output(kwargs['fout'], max_pval)
    elif cltype=='mcl':
        cxs = read_mcl_output(kwargs['fout'])
        pvals,cx_details=None,None
    else:
        assert False, "Wrong cluster type: %s" %cltype
    return Struct(cxs=cxs, pvals=pvals, cx_details=cx_details)

def read_clust_output(fname, max_pval):
    lines = ut.load_list(fname)[1:]
    cxs = [set(line.split(',')[7].strip('"').split()) for line in lines]
    pvals = [line.split(',')[6] for line in lines]
    details = [','.join(line.split(',')[:7]) for line in lines]
    print "Retaining complexes with p < %0.2f." % max_pval
    cxs,pvals,details = keep_pvals(cxs,pvals,details, max_pval)
    return cxs,pvals,details

def keep_pvals(cxs,pvals,details, max_pval):
    inds = [i for i,p in enumerate(pvals) if float(p)<max_pval]
    cxs,pvals,details = [ut.list_inds(lst, inds) for lst in cxs,pvals,details]
    return cxs,pvals,details

def read_mcl_output(fname):
    cxs = [set(c) for c in ut.load_list_of_lists(fname)]
    return cxs

def addpval(cstr,pval):
    cstr = ut.struct_copy(cstr) 
    cstr.params += ',pval_filt=%0.2f'%pval 
    cstr.cxs,cstr.pvals,cstr.cx_details = keep_pvals(cstr.cxs,cstr.pvals,cstr.cx_details,pval) 
    return cstr 

def export_cxs(tested, fname, negmult):
    ut.write_tab_file([(t[0], t[1], ut.rescale(float(t[2]),negmult)) for t in
        tested], fname)

def shell_call(command):
    # http://www.doughellmann.com/PyMOTW/subprocess/
    # use Popen when you want lots of control over the process and output
    #return subprocess.Popen(command,  shell=True, stdout=subprocess.PIPE, )
    # use call when you just want it executed and the output streamed. returns exit code.
    return subprocess.call(command,  shell=True)
    # use check_call when you want an exception generated if it's not successful
    #return subprocess.check_call(command,  shell=True)
    # use check_output to return the output instead of streaming it
    #return subprocess.check_output(command, shell=True)

def merge_maps(cxppis_list, threshold=.5):
    pds = [pd.PairDict(p) for p in cxppis_list]
    pdall = reduce(lambda x,y: pd.pd_combine_ppis(x,y,comb_func=max), pds)
    ps_merged = [list(p)+list(v) for p,v in pdall.d.items() if len([1 for
        whichpd in pds if whichpd.find(p)]) >= len(cxppis_list) * threshold]
    return ps_merged

def remove_dupes(list_cxs, cutoff, func=cp.bader_score, npasses=1):
    """
    cxs: list of sets
    Requires multiple passes to be more certain of fewer complex pairs passing
    the merge threshold. But higher passes results in overcollapse when cutoff
    < 1.
    """
    cxs = list(list_cxs)
    cxs_remaining = list(cxs)
    print "deduping"
    for npass in range(npasses):
        deduped = []
        for i,cx1 in enumerate(cxs):
            #print i
            if cx1 in cxs_remaining:
                cxs_remaining.remove(cx1)
                if (len(cxs_remaining) > 0 and cp.best_match(cx1,
                    cxs_remaining, func=func) >= cutoff):
                    cx2 = cp.best_match_item(cx1,cxs_remaining,func=func)
                    cx_merged = set.union(cx1,cx2)
                    deduped.append(cx_merged)
                    cxs_remaining.remove(cx2)
                    #print "merged"
                else:
                    deduped.append(cx1)
                    #print "not merged"
        cxs = list(deduped)
        cxs_remaining = list(deduped)
    return cxs

def random_complexes(cxs,ps):
    return [set(random.sample(ps, len(c))) for c in cxs]

def random_cxstruct(cxstruct, ps, allppis):
    new_cxst = ut.struct_copy(cxstruct)
    new_cxst.cxs = random_complexes(cxstruct.cxs,ps)
    new_cxst.cxppis = _filter_ints(allppis, new_cxst.cxs)
    return new_cxst

def remove_subcomplexes(cxs, max_size=3):
    """
    Given a list of complexes, return a list without any complexes sized
    max_size or lower that are a subset of another complex.
    """
    cxs_copy = list(cxs)
    return [c for c in cxs if len(c) > max_size or len([c2 for c2 in cxs_copy
        if c.issubset(c2) and not c==c2]) == 0]
