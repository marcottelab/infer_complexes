from __future__ import division
import numpy as np
import subprocess
import os
import complex as co
import utils as ut
import pairdict

def filter_c1(tested, negmult=50, **kwargs):
    cxs = cluster_one(tested, negmult, **kwargs)
    pairs =co.pairs_from_complexes(dict([(i,set(cxs[i]))
                                         for i in range(len(cxs))]))
    pd = pairdict.PairDict(pairs)
    return _filter_ints(tested, pd)
    
def _filter_ints(inlist, pd):
    return [tup for tup in inlist if pd.contains((tup[0],tup[1]))]
    
c1defaults = {
    'min_size': 2,
    'min_density': 0.2,
    'haircut': 0.0,
    'penalty': 2.0,
    'fluff': '--no-fluff',
    'c1path': os.path.expanduser('~/Dropbox/complex/tools')+'/cluster_one.jar',
    'fin': 'temp_interactions.tab',
    'fout': 'temp_complexes.tab'
    }

def cluster_one(tested, negmult, **kwargs):
    kwargs = set_defaults(kwargs, c1defaults)
    export_ints(tested, kwargs['fin'], negmult)
    command = 'java -jar %(c1path)s -s %(min_size)s -d %(min_density)s --haircut %(haircut)s --penalty %(penalty)s %(fluff)s %(fin)s > %(fout)s ' % kwargs
    print command
    shell_call(command)
    ut.printnow('finished cluster_one')
    cxs = ut.load_list_of_lists(kwargs['fout'])
    return cxs

def set_defaults(d, defaultd):
    for k,v in defaultd.items():
        if k not in d:
            d[k] = v
    return d
    
def export_ints(tested, fname, negmult):
    ut.write_tab_file([(t[0], t[1], rescale(t[2],negmult)) for t in
        tested], fname)
    
def rescale(p,n):
    """
    Rescale posterior probability p according to multiple of negatives n.
    """
    return p/(1+(1-p)*n)

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
