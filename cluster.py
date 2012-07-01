from __future__ import division
import numpy as np
import subprocess
import os
import complex as co
import utils as ut
import pairdict

def filter_c1(tested, negmult=50):
    cxs = cluster_one(tested, negmult)
    pairs =co.pairs_from_complexes(dict([(i,set(cxs[i]))
                                         for i in range(len(cxs))]))
    pd = pairdict.PairDict(pairs)
    return _filter_ints(tested, pd)
    
def _filter_ints(inlist, pd):
    return [tup for tup in inlist if pd.contains((tup[0],tup[1]))]
    
def cluster_one(tested, negmult):
    fin='temp_interactions.tab'
    fout = 'temp_complexes.tab'
    export_ints(tested, fin, negmult)
    ut.printnow('running cluster_one')
    #shell_call('java -jar %s/cluster_one.jar -s 3 -d 0.0 --haircut 0.15 \
    #--fluff %s > %s' % (os.path.expanduser('~/Dropbox/complex/tools'), fin,
    #fout))
    shell_call('java -jar %s/cluster_one.jar -s 3 -d 0.0 --haircut 0.15 %s > %s' % (os.path.expanduser('~/Dropbox/complex/tools'), fin, fout))
    ut.printnow('finished cluster_one')
    cxs = ut.load_list_of_lists(fout)
    return cxs

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
