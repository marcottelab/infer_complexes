from __future__ import division
import os
import seqs
import numpy as np
import elution as el
import features as fe
import score as sc
import utils as ut
import Pycluster
from pandas import DataFrame

YSCALE = np.log2(100)

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

def elutions_containing_prots(eluts, sp, query_prots, min_count):
    use_eluts = []
    for e in eluts:
        sp_target = ut.shortname(e.filename)[:2]
        f_inds = set(np.array(np.where(np.max(e.normarr,axis=1) >= min_count)[0])[0])
        if len([p for p in query_prots if p in e.baseid2inds]) > 0:
            if max([ind in f_inds for p in query_prots 
                if p in e.baseid2inds for ind in e.baseid2inds[p]])==True:
                use_eluts.append(e)
    return use_eluts

def save_many_bigprofiles(ind_cxs, unnorm_eluts, fpath, **kwargs):
    for count,(i,cx) in enumerate(ind_cxs):
        fname = os.path.join(fpath, '%s.png' %i)
        tempdir = os.path.join(fpath, '%s.tmp' %i)
        if not os.path.exists(fname):
            if ut.temp_placeholder(tempdir):
                print count, fname
                save_bigprofiles(None, cx, unnorm_eluts, fname, **kwargs)
                os.rmdir(tempdir)


def save_bigprofiles(prots, protids, unnorm_eluts, fname, hires_mult=1, **kwargs):
    import plotting as pl
    nplots = plot_bigprofiles(prots, protids, unnorm_eluts, **kwargs)
    fig = pl.gcf()
    nprots = len(prots) if prots else len(protids)
    fig.set_size_inches(20, 4+(nplots/4)*nprots)
    pl.savefig(fname, bbox_inches='tight', dpi=200*hires_mult)
    pl.clf()

def single_array(gids, unnorm_eluts, sp='Hs', min_count=2,
        remove_multi_base=False, norm_rows=False):
    """
    unnorm_eluts: [el.NormElut(f, sp=sp, norm_cols=False, norm_rows=False) for f in fs]
    """
    import plotting as pl
    use_eluts = elutions_containing_prots(unnorm_eluts, sp, gids, min_count)
    print len(use_eluts), "eluts with proteins"
    ncols = sum([e.normarr.shape[1] for e in use_eluts])
    bigarr = np.zeros((len(gids), ncols))
    startcol = 0
    for e in use_eluts:
        freqarr = ut.normalize_fracs(e.normarr, norm_rows=norm_rows)
        temparr = np.zeros((len(gids), freqarr.shape[1]))
        for i, gid in enumerate(gids):
            if gid in e.baseid2inds:
                inds = list(e.baseid2inds[gid])
                rows = freqarr[inds,:]
                row = np.max(rows, axis=0)
                temparr[i,:] = row
        frac_max = np.max(temparr)
        temparr = np.clip(np.log2(temparr*100 / frac_max), 0, 10)
        bigarr[:, startcol:startcol+freqarr.shape[1]] = temparr
        startcol += freqarr.shape[1]
    return bigarr

def cluster_ids(gids, unnorm_eluts, gt=None, dist='cosine', do_plot=True,
        norm_rows=True, bigarr=None, **kwargs):
    import plotting as pl
    import hcluster
    arr = (bigarr if bigarr is not None else single_array(gids, unnorm_eluts,
        norm_rows=norm_rows))
    ymat = hcluster.pdist(arr, metric=dist)
    zmat = hcluster.linkage(ymat)
    zmat = np.clip(zmat, 0, 10**8)
    if do_plot: pl.figure()
    order = hcluster.dendrogram(zmat, no_plot=bool(1-do_plot), 
            **kwargs)['leaves']
    if do_plot: 
        ax = pl.gca()
        ax.axes.set_xticklabels([gt.id2name[gids[ind]] for ind in order])
        pl.figure() 
        pl.imshow(arr[order,:])
    return list(np.array(list(gids))[order])

def plot_bigprofiles(prots, pids, unnorm_eluts, sp='Hs', min_count=2,
        remove_multi_base=False, gt=None, eluts_per_plot=10,
        do_cluster=True, label_trans=None, do_plot_tree=False, **kwargs):
    """
    supply EITHER prots OR protids, set other to None
    unnorm_eluts: [el.NormElut(f, sp=sp, norm_cols=False, norm_rows=False) for f in fs]
    """
    import plotting as pl
    if prots is not None:
        pids = [gt.name2id[p] for p in prots]
    if do_cluster:
        print "clustering"
        pids = cluster_ids(pids, unnorm_eluts, gt=gt, do_plot=do_plot_tree, 
                **kwargs)
    if gt is not None:
        prots = [gt.id2name[pid] for pid in pids if pid in gt.id2name] #re-order to match
    else:
        prots = pids
        print "No gene names provided--labeling with ids."
    if label_trans: 
        print "Translating names for display."
        # Translate displayed names from base ids according to provided dict
        #prots = [gt.id2name[pid] for pid in pids]
        prots = [label_trans[p] if p in label_trans else p for p in prots]
    prots.reverse(); pids.reverse(); # put them top to bottom
    print "%s proteins" % len(pids)
    use_eluts = elutions_containing_prots(unnorm_eluts, sp, pids, min_count)
    nplots = int(np.ceil(len(use_eluts) / eluts_per_plot))
    maxfracs = 0
    for iplot in range(nplots):
        pl.subplot(nplots, 1, iplot+1)
        plot_eluts = use_eluts[iplot*eluts_per_plot: (iplot+1)*eluts_per_plot]
        startcols = [0]
        for i,e in enumerate(plot_eluts):
            freqarr = ut.normalize_fracs(e.normarr, norm_rows=False)
            sp_target = ut.shortname(e.filename)[:2]
            protsmax = max([np.max(freqarr[r]) for p in pids if p in
                e.baseid2inds for r in e.baseid2inds[p]])
            plot_big_single(freqarr, pids, e.baseid2inds, protsmax,
                    startcols[-1])
            startcols.append(startcols[-1]+freqarr.shape[1])
        label_ys(prots)
        label_xs(startcols, [ut.shortname(e.filename) for e in plot_eluts])
        pl.grid(False)
        maxfracs = maxfracs if maxfracs > startcols[-1] else startcols[-1]
    for iplot in range(nplots):
        pl.subplot(nplots, 1, iplot+1)
        pl.xlim(0,maxfracs)
    pl.subplots_adjust(hspace=5/len(prots))
    return nplots

def label_xs(lefts, labels):
    import plotting as pl
    for left in lefts:
        pl.axvline(x=left, linewidth=.5, color='gray')
    ax = pl.gca()
    ax.axes.set_xticks([x+2 for x in lefts])
    ax.axes.set_xticklabels(labels)
    labels = ax.get_xticklabels() 
    for label in labels: 
        label.set_rotation(30) 

def label_ys(labels):
    import plotting as pl
    ax = pl.gca()
    bottoms = [YSCALE*i for i in range(len(labels))]
    for b in bottoms: pl.axhline(y=b, linewidth=.5, color='gray')
    ax.axes.set_yticks([b+YSCALE/2 for b in bottoms])
    ax.axes.set_yticklabels(labels)

def plot_big_single(arr, pids, baseid2inds, maxcount, startcol):
    import plotting as pl
    for i,pid in enumerate(pids):
        if pid in baseid2inds:
            for rowid in baseid2inds[pid]:
                row = arr[rowid]
                bottom = YSCALE*i
                plot_vals = np.clip(np.log2([x*100/maxcount for x in
                    row]),0,100)
                pl.bar([x+startcol for x in range(len(row))], plot_vals,
                        color=pl.COLORS[i%len(pl.COLORS)],
                        align='center',width=1,linewidth=0, bottom=bottom,
                        antialiased=False)

def plot_profiles(prots, eluts, sp='Hs', plot_sums=True, shape=None,
        min_count=2):
    """
    shape: (m,n) = m rows, n columns
    eluts: [el.NormElut(f, sp, norm_rows=False, norm_cols=False) for f in
    fs]
    """
    import plotting as pl
    gt = seqs.GTrans()
    use_eluts = elutions_containing_prots(eluts, sp, seqs.names2ids(prots),
            min_count)
    shape = shape if shape else ut.sqrt_shape(len(use_eluts)+1)
    fig = pl.figure()
    for i,e in enumerate(use_eluts):
        sp_target = ut.shortname(e.filename)[:2]
        pl.subplot(shape[0],shape[1],i+1)
        pl.title(ut.shortname(e.filename))
        pids = [gt.name2id[p] for p in prots]
        protsmax = max([np.max(e.normarr[r]) for p in pids if p in e.baseid2inds for
            r in e.baseid2inds[p]])
        plot_prots(e, pids, e.baseid2inds, protsmax)
        if plot_sums:
            # plot total spectral counts normalized to match biggest peak
            sums = np.sum(e.normarr,axis=0)
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
                row = elut.normarr[rowid]
                bottom = np.log2(maxcount)*i
                pl.bar(range(row.shape[1]),
                        np.clip(np.log2(row[0,:].T+.1),0,1000),
                        color=pl.COLORS[i%len(pl.COLORS)], align='center',width=1,linewidth=0,
                        bottom=bottom)
    pl.xlim(0,row.shape[1])
    pl.ylim(-.1, np.log2(maxcount)*len(pids))

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

def treeview_eluts(name, fs, base_sp='Hs', gt=None, ids=None, bigarr=None):
    print "Building big array."
    ids = ids or el.all_prots(fs,sp_base=base_sp)
    gt = gt or seqs.GTrans()
    gene_names = [gt.id2name.get(x,x) for x in ids]
    if bigarr is None:
        unes = [el.NormElut(f,sp_base=base_sp, norm_rows=False,norm_cols=False)
                for f in fs]
        bigarr = single_array(ids, unes, sp=base_sp, norm_rows=True)
        cols = ut.flatten([[ut.shortname(u.filename)+'_'+str(i) 
            for i in range(u.normarr.shape[1])] for u in unes])
        df = DataFrame(data=bigarr, index=list(gene_names), columns=cols)
    print "Writing to temp file and reading into Pycluster object."
    tempf = name + '.tab'
    df.to_csv(tempf, sep='\t')
    record = Pycluster.Record(open(tempf))
    os.remove(tempf)
    print "Clustering."
    tree = record.treecluster()
    print "Saving."
    record.save(name, tree)

def gene_ppis(gnames, gids, ppis, gtrans=None, sp='Hs'):
    """
    Provide pres.ppis as input, get all interactors with the provied genes via
    gnames or gids.
    """
    gt = gtrans or seqs.GTrans(sp=sp)
    gids = gids or [gt.name2id[n] for n in gnames]
    gids = set(gids)
    return [(i,(gt.id2name[p[0]], gt.id2name[p[1]], p[2], gt.id2desc[p[0]],
        gt.id2desc[p[1]])) for i,p in enumerate(ppis) if p[0] in gids or p[1] in gids]

def ppis_fracspassing_counts(ppis, obs, exclude_ppis=None, cutoff=0.5):
    """
    For a limited set of ppis (say top 15k), return a list that also includes
    as the final column the number of fractionations in which that ppi passes
    the threshold.
    exclude_ppis is chiefly for excluding cxppis.
    ppis should have been generated from obs, so pair order should be the same,
    although list order is different.
    """
    def pair2ind(items):
        return ut.list_inv_to_dict(((x[0],x[1]) for x in items))
    dobs = pair2ind(obs)
    obs_forppis = obs[[dobs[(p[0],p[1])] for p in ppis]]
    npassing_obs_forppis = fe.passing_fractionations(obs_forppis)
    newppis = []
    for i,p in enumerate(ppis):
        newppis.append(tuple(list(p) + [npassing_obs_forppis[i]])) 
    if exclude_ppis is not None:
        d_exclude = pair2ind(exclude_ppis)
        newppis = [p for p in newppis if (p[0],p[1]) not in d_exclude]
    print "NOT SORTED!"
    return newppis

