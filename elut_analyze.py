from __future__ import division
from pylab import *
import elution as el
import utils as ut
import orth
import collections
import scipy

def prot_counts(fs, min_count=2):
    """
    Sum up all the spectral counts for all the proteins in a set of
    fractionations.  
    Filtered s.t. any returned protein will have at least min_count counts in
    one fraction of one of the fractionations.
    Return a dict: {prot1:count1, prot2:count2, ...}
    """
    allprots = el.all_prots(fs, min_count=min_count)
    pcounts = collections.defaultdict(int)
    for f in fs:
        e = el.load_elution(f)
        psums = np.sum(np.array(e.mat),axis=1)
        for p,psum in zip(e.prots,psums):
            if p in allprots:
                pcounts[p] += psum
    return pcounts

def prot_conservation(fs,sp1,sp2, gridsize=30, **kwargs):
    """
    Currently only uses 1 to 1 orthologs, so odict should be a simple flat dict
    of genesa:genesb.
    """
    fs1,fs2 = [[f for f in fs if ut.shortname(f)[:2]==sp] for sp in sp1,sp2]
    odict = orth.odict_1to1(sp1,sp2)
    pcounts = [prot_counts(fs) for fs in fs1,fs2]
    pc1,pc2 = zip(*[(pcounts[0][p], pcounts[1][odict[p]]) 
        for p in pcounts[0] if p in odict])
    logpc1,logpc2 = [np.log2(pc) for pc in pc1,pc2]
    hexbin(logpc1,logpc2,gridsize=gridsize,**kwargs)
    xlabel('%s log2 unique spectral counts' %sp1)
    ylabel('%s log2 unique spectral counts' %sp2)
    title('%s-%s: R^2 = %0.2f, %s 1-1 ortholog pairs with %s data' % (sp1,sp2,
        scipy.stats.pearsonr(pc1,pc2)[0]**2, len(pc1), sp1))
    return pc1,pc2

def plot_conservations(fs,sp,othersps,**kwargs):
    for i,osp in enumerate(othersps):
        subplot(2,np.ceil(len(othersps)/2),i+1)
        prot_conservation(fs, sp, osp, **kwargs)

def hist_prot_counts(fs, **kwargs):
    pcounts = ut.i1(prot_counts(fs).items())
    hist(np.log2(pcounts), **kwargs)


