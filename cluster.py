from __future__ import division
import numpy as np
import subprocess
import os
import corum as co
import utils as ut
import pairdict
import operator

mcl_command = '%(mclpath)s %(fin)s -I %(I)s -o %(fout)s --abc'
c1_command = 'java -jar %(c1path)s -s %(min_size)s -d %(min_density)s --haircut %(haircut)s --penalty %(penalty)s --max-overlap %(max_overlap)s %(fluff)s %(fin)s > %(fout)s '

def multi_clust(tested, fracs=[0.03,.05,.07], d_lgm=[.03,.05,.07],
        m_lgm=[40,60,100], d_smm=[.2,.3,.4,.45,.5], m_smm=[5,10,20,40],
        **kwargs):
    clusts = reduce(operator.add,
            [[('density=%s,frac=%s,negmult=%s'%(d,f,m),
                filter_clust(tested[:int(f*len(tested))], negmult=m,
                    min_density=d))
                for d in ds for f in fs for m in ms]
                for ds,fs,ms in (d_lgm,fracs,m_lgm),(d_smm,fracs,m_smm)])
    return clusts

def filter_clust(tested, negmult=50, command=c1_command, do_dedupe=True,
        **kwargs):
    cxs = cluster(tested, negmult, command, **kwargs)
    if do_dedupe: cxs = remove_dupes(cxs)
    return cxs, _filter_ints(tested, cxs)
    
def _filter_ints(inlist, cxs):
    d_cxs = dict([(i,set(c)) for i,c in enumerate(cxs)]) 
    pairs = co.pairs_from_complexes(d_cxs)
    pd = pairdict.PairDict(pairs)
    return [tup for tup in inlist if pd.contains((tup[0],tup[1]))]

clust_defaults = {
    'I': 1.8,
    'min_size': 2,
    'min_density': 0.2,
    'haircut': 0.0,
    'penalty': 2.0,
    'max_overlap': 0.6,
    'fluff': '--no-fluff',
    'c1path': os.path.expanduser('~/Dropbox/complex/tools')+'/cluster_one.jar',
    'mclpath': os.path.expanduser('~/local/bin')+'/mcl',
    'fin':
        os.path.expanduser('~/bigdata/complex/predictions')+
        '/temp_interactions.tab',
    'fout': os.path.expanduser('~/bigdata/complex/predictions') +
        '/temp_complexes.tab'
    }

def thresh(tested, score):
    return tested[:n_thresh(tested, score)]

def n_thresh(tested, score):
    for i,t in enumerate(tested):
        if t[2] < score:
            return i
    
def cluster(tested, negmult, command, **kwargs):
    kwargs = set_defaults(kwargs, clust_defaults)
    export_cxs(tested, kwargs['fin'], negmult)
    print command
    shell_call(command % kwargs)
    ut.printnow('finished clustering')
    cxs = [set(c) for c in ut.load_list_of_lists(kwargs['fout'])]
    return cxs

def set_defaults(d, defaultd):
    for k,v in defaultd.items():
        if k not in d:
            d[k] = v
    return d
    
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


def remove_dupes(cxs):
    """
    cxs: list of sets
    """
    deduped = []
    for c in cxs:
        if c not in deduped:
            deduped.append(c)
    return deduped
