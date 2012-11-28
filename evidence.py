from __future__ import division
import utils as ut
import seqs
import numpy as np
import elution as el
import score as sc

def all_ev(gnames, arr_ev, scores, req_gnames=set([])):
    """
    Given a list of gene names, return the rows of evidence for any pairs.
    Filter only to those pairs present in scores (weighted predicted
    interaction list--probably the cx_filtered version), if provided.
    """
    gt = seqs.GTrans()
    ids = set([gt.name2id[n] for n in gnames])
    if scores:
        # don't have to check flip(s[0],s[1]) since scores derived from arr_ev
        pairs = dict([((s[0],s[1]),s[2]) for s in scores 
            if s[0] in ids and s[1] in ids])
        allevs, titles = ([tuple([gt.id2name[i] for i in r[0],r[1]] +
            [pairs[(r[0],r[1])]] + list(r)[3:]) for r in arr_ev if (r[0],r[1])
            in pairs]), ['total_score']+list(arr_ev.dtype.names[3:])
    else: 
        allevs, titles = ([tuple([gt.id2name[i] for i in
            r[0],r[1]]+list(r)[3:]) for r in arr_ev if r[0] in ids and r[1] in
            ids], list(arr_ev.dtype.names[3:]))
    if req_gnames:
        allevs = [a for a in allevs if a[0] in req_gnames or
            a[1] in req_gnames]
    return allevs, titles

def ev_output(gnames, arr_ev, scores=None, **kwargs):
    ev_lot, titles = all_ev(gnames, arr_ev, scores, **kwargs)
    labels = ['gene1','gene2'] + titles
    return zip(labels, *ev_lot)

def disp_ev(gnames, arr_ev, scores=None, cutoff=.25, dispn=10, **kwargs):
    evs = ev_output(gnames, arr_ev, scores=scores, **kwargs)
    headers = evs[:2]
    data = evs[2:]
    data.sort(key=lambda(x):len([v for v in x[1:] if v>cutoff]), reverse=True)
    output = headers+data
    for r in output[:dispn]: print r
    return evs

def elutions_containing_prots(fs, sp, query_prots, min_count,
        remove_multi_base):
    usefs = []
    for f in fs:
        e = el.load_elution(f)
        sp_target = ut.shortname(f)[:2]
        baseid2inds = sc.orth_indices(sp, sp_target, e.prots, remove_multi_base)
        f_inds = set(np.array(np.where(np.max(e.mat,axis=1) >= min_count)[0])[0])
        if len([p for p in query_prots if p in baseid2inds]) > 0:
            if max([ind in f_inds for p in query_prots 
                if p in baseid2inds for ind in baseid2inds[p]])==True:
                usefs.append(f)
    return usefs

def plot_profiles(prots, fs, sp='Hs', plot_sums=True, shape=None, min_count=2,
        remove_multi_base=False):
    """
    shape: (m,n) = m rows, n columns
    """
    import plotting as pl
    gt = seqs.GTrans()
    usefs = elutions_containing_prots(fs, sp, seqs.names2ids(prots), min_count,
            remove_multi_base)
    shape = shape if shape else ut.sqrt_shape(len(usefs)+1)
    fig = pl.figure()
    for i,f in enumerate(usefs):
        e = el.load_elution(f)
        sp_target = ut.shortname(f)[:2]
        baseid2inds = sc.orth_indices(sp, sp_target, e.prots, remove_multi_base)
        pl.subplot(shape[0],shape[1],i+1)
        pl.title(ut.shortname(f))
        pids = [gt.name2id[p] for p in prots]
        protsmax = max([np.max(e.mat[r]) for p in pids if p in baseid2inds for
            r in baseid2inds[p]])
        plot_prots(e, pids, baseid2inds, protsmax)
        if plot_sums:
            # plot total spectral counts normalized to match biggest peak
            sums = np.sum(e.mat,axis=0)
            fmax = np.max(sums)
            pl.plot(range(sums.shape[1]),
                    np.log2(sums[0,:]).T*np.log2(protsmax)*len(pids)/np.log2(fmax), 
                    color='k', linestyle='-', linewidth=.5)
    # make legend with all prots
    pl.subplot(shape[0],shape[1],0)
    for p in prots: pl.plot(0,label=p)
    pl.legend()

def plot_prots(elut, pids, baseid2inds, maxcount):
    import plotting as pl
    for i,pid in enumerate(pids):
        if pid in baseid2inds:
            for rowid in baseid2inds[pid]:
                row = elut.mat[rowid]
                bottom = np.log2(maxcount)*i
                pl.bar(range(row.shape[1]),
                        np.clip(np.log2(row[0,:].T+.1),0,1000),
                        color=pl.COLORS[i], align='center',width=1,linewidth=0,
                        bottom=bottom)
                pl.xlim(0,row.shape[1])

#def plot_profiles(prots, fs, sp='Hs', plot_sums=True, shape=None, min_count=2,
        #remove_multi_base=False):
    #"""
    #shape: (m,n) = m rows, n columns
    #"""
    #import plotting as pl
    #gt = seqs.GTrans()
    #usefs = elutions_containing_prots(fs, sp, seqs.names2ids(prots), min_count,
            #remove_multi_base)
    #shape = shape if shape else ut.sqrt_shape(len(usefs)+1)
    #fig = pl.figure()
    #for i,f in enumerate(usefs):
        #e = el.load_elution(f)
        #sp_target = ut.shortname(f)[:2]
        #baseid2inds = sc.orth_indices(sp, sp_target, e.prots, remove_multi_base)
        #ax = fig.add_subplot(shape[0],shape[1],i+1)
        #pl.title(ut.shortname(f))
        ##protsmax = 0
        #for i,(pid,pname) in enumerate([(gt.name2id[p],p) for p in prots]):
            #if pid in baseid2inds:
                #plot_prot_separately(baseid2inds[pid], e, fig, ax, len(prots)+1, i) 
                ##plot_inds_separately(baseid2inds[pid], e, ax)
        ## plot total spectral counts normalized to match biggest peak
        #sums = np.sum(e.mat,axis=0)
        #fmax = np.max(sums)
        ##pl.plot(range(sums.shape[1]), sums[0,:].T*protsmax/fmax, color='k',
                ##linestyle=':')
        #newax = fig.add_axes(hor_slice_position(ax.get_position(),
            #len(prots)+1, len(prots)))
        #newax.plot(range(sums.shape[1]), sums[0,:].T, color='k', linestyle=':')
    ## make legend with all prots
    #pl.subplot(shape[0],shape[1],0)
    #for p in prots: pl.plot(0,label=p)
    #pl.legend()

def plot_inds_together(inds, elut):
    for rowid in baseid2inds[pid]:
        row = elut.mat[rowid]
        if np.max(row) >= min_count:
            rowmax = np.max(row)
            protsmax = rowmax if rowmax > protsmax else protsmax
            pl.plot(range(row.shape[1]),row[0,:].T, label=pname,
                    color=pl.COLORS[i],linewidth=1)

def plot_prot_separately(inds, elut, fig, ax, nprots, index):
    import plotting as pl
    newax = fig.add_axes(hor_slice_position(ax.get_position(), nprots, index))
    for ind in inds:
        row = elut.mat[ind]
        newax.plot(range(row.shape[1]),row[0,:].T, color=pl.COLORS[index],
                linewidth=2)
    xlim(0,row.shape[1])

def hor_slice_position(bbox, nslices, index):
    """
    new rect position: left, bottom, width, height
    """
    newheight = bbox.height / nslices
    newymin = bbox.ymin + index * newheight
    newbox = (bbox.xmin, newymin, bbox.width, newheight)
    return newbox

def plot_sums(fs, shape=None):
    import plotting as pl
    shape = shape if shape else ut.sqrt_shape(len(fs))
    for i,f in enumerate(fs):
        e = el.load_elution(f)
        pl.subplot(shape[0],shape[1],i+1)
        pl.title(ut.shortname(f))
        sums = np.sum(e.mat,axis=0)
        pl.plot(range(sums.shape[1]), sums[0,:].T)
