from __future__ import division
import numpy as np
import subprocess
import os
import corum as co
import utils as ut
import pairdict as pd
from Struct import Struct
import operator
import random

mcl_command = '%(mclpath)s %(fin)s -I %(I)s -o %(fout)s --abc'
c1_command = 'java -jar %(c1path)s -s %(min_size)s -d %(min_density)s --haircut %(haircut)s --penalty %(penalty)s --seed-method edges -F csv --max-overlap %(max_overlap)s %(fluff)s %(fin)s > %(fout)s '


def multi_clust(tested, fracs=[.01,.02,.03],
        d_sets=[[.02,.05,.07],[.003,.006, .01,.03]],
        m_sets=[[1,3,5,10],[10,15,25,50]], penalties=[1], overlaps=[.7],
        max_pval=1, savef=None, runid=random.randrange(0,100), **kwargs):
    print "random id for these runs:", runid
    clusts = []
    nsets = len(d_sets)
    f_sets = [fracs for i in range(nsets)]
    p_sets = [penalties for i in range(nsets)]
    o_sets = [overlaps for i in range(nsets)]
    for fs,ds,ms,ps,ovs in zip(f_sets, d_sets, m_sets, p_sets, o_sets):
        for f in fs:
            for d in ds:
                for m in ms:
                    for p in ps:
                        for o in ovs:
                            cxstruct = filter_clust(tested[:int(f*len(tested))],
                                    negmult=m, min_density=d, runid=runid,
                                    penalty=p, max_pval=max_pval,
                                    max_overlap=o, **kwargs)
                            cxstruct.params = ('density=%s,frac=%s,negmult=%s,penalty=%s,max_pval=%s,max_overlap=%s' % (d,f,m,p, max_pval, o))
                            clusts.append(cxstruct)
                            if savef and (len(clusts) % 10 == 1):
                                ut.savepy(clusts, ut.pre_ext(savef, "clusts_%s" % runid))
    #clusts = reduce(operator.add,
    return clusts, runid

def filter_clust(tested, negmult=50, command=c1_command, do_dedupe=True,
        **kwargs):
    print "Negatives rescale factor:", negmult
    cxstruct = cluster(tested, negmult, command, **kwargs)
    if do_dedupe: cxstruct.cxs = remove_dupes(cxstruct.cxs)
    cxstruct.cxppis = _filter_ints(tested, cxstruct.cxs)
    print "%s Deduped complexes; %s cxppis" % (len(cxstruct.cxs),
            len(cxstruct.cxppis))
    return cxstruct
    
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
    
def cluster(tested, negmult, command, max_pval=0.2, **kwargs):
    kwargs = ut.dict_set_defaults(kwargs, clust_defaults)
    if 'runid' in kwargs: # keep temp files separate
        runid = str(kwargs['runid']) 
        kwargs['fin'] = ut.pre_ext(kwargs['fin'], runid)
        kwargs['fout'] = ut.pre_ext(kwargs['fout'], runid)
    export_cxs(tested, kwargs['fin'],negmult)
    command = command % kwargs
    print command
    shell_call(command)
    cxs, pvals, cx_details = read_clust_output(kwargs['fout'], max_pval)
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
    keep = [(c,p,d) for c,p,d in zip(cxs,pvals,details) if float(p)<max_pval]
    cxs,pvals,details = zip(*keep)
    return cxs,pvals,details

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

def remove_dupes(cxs):
    """
    cxs: list of sets
    """
    deduped = []
    for c in cxs:
        if c not in deduped:
            deduped.append(c)
    return deduped
