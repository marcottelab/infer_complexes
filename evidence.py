from __future__ import division
import utils as ut
import seqs
import numpy as np
import elution as el
import score as sc

def all_ev(gnames, arr_ev, scores):
    """
    Given a list of gene names, return the rows of evidence for any pairs.
    Filter only to those pairs present in scores (weighted predicted
    interaction list--probably the cx_filtered version), if provided.
    """
    gt = seqs.GTrans()
    ids = set([gt.name2id[n] for n in gnames])
    if scores:
        pairs = dict([((s[0],s[1]),s[2]) for s in scores 
            if s[0] in ids and s[1] in ids])
        return ([tuple([gt.id2name[i] for i in r[0],r[1]] +
            [pairs[(r[0],r[1])]] + list(r)[3:]) for r in arr_ev if (r[0],r[1])
            in pairs]), ['total_score']+list(arr_ev.dtype.names[3:])
    else: 
        return ([tuple([gt.id2name[i] for i in r[0],r[1]]+list(r)[3:]) for r in
                arr_ev if r[0] in ids and r[1] in ids],
                list(arr_ev.dtype.names[3:]))

def ev_output(gnames, arr_ev, scores=None):
    ev_lot, titles = all_ev(gnames, arr_ev, scores)
    labels = ['gene1','gene2'] + titles
    return zip(labels, *ev_lot)

def disp_ev(gnames, arr_ev, scores=None, dispn=10):
    evs = ev_output(gnames, arr_ev, scores=scores)
    headers = evs[:2]
    data = evs[2:]
    data.sort(key=lambda(x):len([v for v in x[1:] if v>.5]), reverse=True)
    output = headers+data
    for r in output[:dispn]: print r
    return output

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
    for i,f in enumerate(usefs):
        e = el.load_elution(f)
        sp_target = ut.shortname(f)[:2]
        baseid2inds = sc.orth_indices(sp, sp_target, e.prots, remove_multi_base)
        pl.subplot(shape[0],shape[1],i+1)
        pl.title(ut.shortname(f))
        protsmax = 0
        for i,(pid,pname) in enumerate([(gt.name2id[p],p) for p in prots]):
            if pid in baseid2inds:
                for rowid in baseid2inds[pid]:
                    row = e.mat[rowid]
                    if np.max(row) >= min_count:
                        rowmax = np.max(row)
                        protsmax = rowmax if rowmax > protsmax else protsmax
                        pl.plot(range(row.shape[1]),row[0,:].T, label=pname,
                                color=pl.COLORS[i],linewidth=1)
        # plot total spectral counts normalized to match biggest peak
        sums = np.sum(e.mat,axis=0)
        fmax = np.max(sums)
        pl.plot(range(sums.shape[1]), sums[0,:].T*protsmax/fmax, color='k',
                linestyle=':')
    # make legend with all prots
    pl.subplot(shape[0],shape[1],0)
    for p in prots: pl.plot(0,label=p)
    pl.legend()


def plot_sums(fs, shape=None):
    import plotting as pl
    shape = shape if shape else ut.sqrt_shape(len(fs))
    for i,f in enumerate(fs):
        e = el.load_elution(f)
        pl.subplot(shape[0],shape[1],i+1)
        pl.title(ut.shortname(f))
        sums = np.sum(e.mat,axis=0)
        pl.plot(range(sums.shape[1]), sums[0,:].T)
