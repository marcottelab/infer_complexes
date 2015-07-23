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
    pcounts = collections.defaultdict(float)
    for f in fs:
        e = el.load_elution(f)
        psums = np.sum(np.array(e.mat),axis=1)
        frac_sum = sum(psums)
        norm_term = 1 / (frac_sum * len(fs))
        for p,psum in zip(e.prots,psums):
            if p in allprots:
                pcounts[p] += (psum * norm_term)
    return pcounts


def prot_conservation(fs,sp1,sp2, gridsize=30, od11=None, return_data=False,
        filter_both=True, use_title=True, extent=[-22,-6,-22,-6], fontsize=18,
        **kwargs):
    """
    Currently only uses 1 to 1 orthologs, so odict should be a simple flat dict
    of genesa:genesb.
    """
    if sp1==sp2:
        return
    fs1,fs2 = [[f for f in fs if ut.shortname(f)[:2]==sp] for sp in sp1,sp2]
    odict = orth.odict_1to1(sp1,sp2) if od11 == None else od11
    pc1_all,pc2_all = [prot_counts(fs) for fs in fs1,fs2]
    if filter_both:
        ps_use = [p for p in odict if (pc1_all[p]>0 and pc2_all[odict[p]]>0)]
    else:
        ps_use = [p for p in pc1_all if p in odict]
    pc1,pc2 = zip(*[(pc1_all[p], pc2_all[odict[p]]) for p in ps_use])
    logpc1,logpc2 = [np.log2(pc) for pc in pc1,pc2]
    plot(extent[:2],extent[2:],'k:', linewidth=1)
    hexbin(logpc1,logpc2,gridsize=gridsize,**kwargs)
    #if use_title:
        #xlabel('%s log2 unique spectral counts' %sp1)
        #ylabel('%s log2 unique spectral counts' %sp2)
        #title('%s-%s: spearmanR: %0.2f, %s 1-1 nonzero ortholog pairs' %
                #(sp1,sp2, scipy.stats.spearmanr(pc1,pc2)[0], len(pc1)))
    rval = scipy.stats.spearmanr(pc1,pc2)[0]
    annotate('R=%0.2f\nN=%s' % (rval, len(pc1)), xy=(.05,.7),
            xycoords=("axes fraction"), fontsize=fontsize)
    if return_data:
        return pc1,pc2

def plot_conservations(fs,sp,othersps,**kwargs):
    for i,osp in enumerate(othersps):
        subplot(2,np.ceil(len(othersps)/2),i+1)
        prot_conservation(fs, sp, osp, **kwargs)
        grid("off")
        gca().xaxis.set_ticklabels([])
        gca().yaxis.set_ticklabels([])
        draw()

def all_by_all_species(fs, sps, extent=[-22,-6,-22,-6], fontsize=18, **kwargs):
    nsps = len(sps)
    kwargs["extent"] = extent
    for j,s1 in enumerate(sps):
        for i,s2 in enumerate(sps):
            print s1, s2, j*nsps+i+1
            subplot(nsps, nsps, j*nsps+i+1)
            if s1!=s2:
                plot(extent[:2],extent[2:],'k:', linewidth=1)
            result = prot_conservation(fs, s1, s2, use_title=False,
                    return_data=True, **kwargs)
            if result:
                pc1,pc2=result
                rval = scipy.stats.spearmanr(pc1,pc2)[0]
                annotate('R=%0.2f\nN=%s' % (rval, len(pc1)), xy=(.05,.7),
                    xycoords=("axes fraction"), fontsize=fontsize)
            grid("off")
            gca().xaxis.set_ticklabels([])
            gca().yaxis.set_ticklabels([])
            draw()
            

def hist_prot_counts(fs, in_prots_sets=[None],**kwargs):
    """
    Usually set linewidth=3, histtype='step', range=(0,18), bins=36
    """
    pcounts = prot_counts(fs).items()
    for prots in in_prots_sets:
        if prots:
            pcounts = [(p,c) for p,c in pcounts if p in prots]
        hist(np.log2(ut.i1(pcounts)), **kwargs)


