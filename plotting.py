from __future__ import division
import numpy as np
import scipy
import pylab
from pylab import *
import random
import itertools
import cv
from Struct import Struct
import utils as ut
import pairdict as pd
import pandas as pa
import ppi_utils as pu


#COLORS = ['#4571A8', 'black', '#A8423F', '#89A64E', '#6E548D', '#3D96AE',
           #'#DB843D', '#91C4D5', '#CE8E8D', '#B6CA93', '#8EA5CB', 'yellow',
           #'gray', 'blue']
#COLORS_BLACK = ['#4571A8', 'white', '#A8423F', '#89A64E', '#6E548D', '#3D96AE',
           #'#DB843D', '#91C4D5', '#CE8E8D', '#B6CA93', '#8EA5CB', 'yellow',
           #'gray', 'blue']
COLORSTRING = "4571A8, 000000, A8423F, 89A64E, 6E548D, 3D96AE, DB843D, 91C4D5, CE8E8D, B6CA93, 8EA5CB, FFFF00, 404040, 0000FF"
COLORS_WHITE = ["#"+c for c in COLORSTRING.split(', ')]
COLORSTRING_BLACK = "4571A8, FFFFFF, A8423F, 89A64E, 6E548D, 3D96AE, DB843D, 91C4D5, CE8E8D, B6CA93, 8EA5CB, FFFF00, 404040, 0000FF"
COLORS_BLACK = ["#"+c for c in COLORSTRING_BLACK.split(', ')]
COLORS = COLORS_WHITE
          

def plot_result(result, ppis=None, **kwargs):
    ppis = ppis if ppis else result.ppis
    kwargs['label'] = kwargs.get('label','') +' '+ result.name
    pr_plot(ppis, result.ntest_pos, **kwargs)

def boot_resample(extr_exte):
    return [Struct(names=ex.names,examples=ut.sample_wr(ex.examples, len(ex.examples))) for ex in extr_exte]

def rolling_scores(tested, true_ints=None, show=1000, window=50, rescale=0,
        **kwargs):
    if rescale > 0:
        tested = [(t[0],t[1],ut.rescale(t[2],rescale), t[3]) for t in tested]
    if true_ints:
        tested = cv.tested_from_trues(tested, true_ints)
    padded = list(np.zeros((50,4)))+list(tested)
    rolling = [len([t for t in padded[i:i+window] if t[3]==1])/window for i in range(show)]
    plot(rolling, **kwargs)
    plot([t[2] for t in tested[:show]], **kwargs)
    xlabel('starting index in scored examples')
    ylabel('fraction true in index:index+%s'%window)
    legend(['fraction true','score'])


def cumulative_precision(tested, true_ints=None, return_data=False, 
        **kwargs):
    if true_ints:
        tested = cv.tested_from_trues(tested, true_ints)
    hits = [t[3] for t in tested]
    precision = [sum(hits[:i+1])/(i+1) for i in range(len(hits))]
    plot(precision, **kwargs)
    plot([t[2] for t in tested], **kwargs)
    xlabel('Index of PPI')
    ylabel('Cumulative Precision; PPI score')
    legend(['Cumulative Precision','PPI score'])
    if return_data: 
        return precision

def bar_plot(names, yvals, slanted_names=False, ax=None, **kwargs):
    fig = figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(names)), yvals, align='center', **kwargs)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    if slanted_names: fig.autofmt_xdate()
    show()
    return ax

def stacked_bar(names, values):
    """
    values is a lol. values[0] corresponds to the values for names[0].
    """
    valuesT = zip(*values)
    padded = [[0]*len(valuesT[0])] + valuesT # 0s in first row is helpful
    arr = np.array(padded)
    arrcum = np.cumsum(arr, axis=0)
    for i in range(1, arr.shape[1]):
        bar(range(1,arr.shape[0]+1), arr[i], align='center', bottom=arrcum[i-1],
                color=COLORS[i-1])
    ax = gca()
    ax.set_xticklabels(['']+names)
    df = pa.DataFrame(arr[1:][::-1], columns=names)
    print df
    return df
    
# def cluster(corr):
#     # corr: a matrix of similarity scores, such as a covariance matrix
#     ymat = hcluster.pdist(corr)
#     zmat = hcluster.linkage(ymat)
#     figure()
#     order = hcluster.dendrogram(zmat)['leaves']
#     figure()
#     imshow(corr[order,:][:,order])
#     # check for failure signs
#     for i in random.sample(range(len(order)),10):
#         if order[i] - order[i-1] == 1:
#             print 'HEY!! probable clustering failure.'
#             break
#     return order

def roc_plot(cvpairs, **kwargs):
    xs,ys = cv.roc(cvpairs) 
    auroc = cv.auroc(xs,ys)
    kwargs['label'] = kwargs.get('label','') + ' %.3f' % auroc
    plot(xs, ys, **kwargs)
    plot([0,xs[-1]], [0,ys[-1]], 'k--')

def pr_plot(cv_pairs, total_trues, rescale=None, style=None, prec_test=None,
        true_ints=None, return_data=False, do_plot=True, **kwargs):
    """
    rescale: adjust precision values assuming rescale times as many negatives
    total_trues:
    - None for just displaying recall count instead of fraction
    - 'auto' to calculate from the supplied tested cv_pairs
    - integer to use that supplied integer as total trues
    """
    if true_ints:
        pdtrues = pd.PairDict(true_ints)
        cv_pairs = [(p[0],p[1],p[2],1 if pdtrues.contains(tuple(p[:2])) else 0) for p
                in cv_pairs]
    if total_trues == 'auto':
        total_trues = len([t for t in cv_pairs if t[3]==1])
    recall,precision = cv.pr(cv_pairs) 
    if rescale:
        precision = [ p / (p + (1-p) * rescale) for p in precision]
    if prec_test:
        kwargs['label'] = kwargs.get('label','') + (' Re:%0.2f' %
        cv.calc_recall(precision,prec_test, total_trues)) + (' @ Pr:%0.2f'
            % prec_test)
    if total_trues:
        recall = [r/total_trues for r in recall]
    args = [style] if style is not None else []
    if do_plot:
        plot(recall, precision, *args, **kwargs)
        xlabel('Recall: TP/(TP+FN)')
        ylabel('Precision: TP/(TP+FP)')
        ylim(-0.02,1.02)
        xlim(xmin=-0.002)
        legend()
    if return_data:
        return recall,precision

def imshow2(*args):
    imshow(*args, interpolation='nearest', aspect='auto',
           cmap='bone', vmin=0)
           # Need vmin = 0 so the lowest values aren't represented 
           # as say black by the colorbar if they're not 0.
           # Colormaps: bone, gray

def examples_dist_arr(arr, score_indices, ncols=1, use_legend=False,
        return_data=False, use_randoms=False,  **kwargs):
    extra = 1 if use_legend else 0
    nplots = len(score_indices)+extra
    pos,neg = pos_neg_from_arr(arr, use_all=use_randoms)
    results = []
    for i, (ind, name) in enumerate([(ind,arr.dtype.names[ind]) for ind in score_indices]):
        subplot(int(np.ceil(nplots/ncols))+extra, ncols, i+1+extra)
        hp,hn = examples_dist_single(pos_neg=(pos,neg), name=name, **kwargs)
        #kwargs['range'] = kwargs['range'] if 'range' in kwargs else \ [func([func(data) for data in [pos,neg]]) for func in [min,max]]
        title(name)
        results.append((hp,hn))
    if use_legend and ((not 'odds' in kwargs) or (not kwargs['odds'])): 
        subplot(int(nplots/ncols)+2,ncols,1)
        legend([hp[0],hn[0]],['Pos','Neg'])
    if return_data:
        return results

def pos_neg_from_arr(arr, use_all=False):
    pos,neg = [arr[[i for i in range(len(arr)) if arr[i][2]==t]] for t in 1,0]
    if use_all:
        print "Using all for negatives."
        neg = arr
    return pos,neg

def examples_dist_single(arr=None, pos_neg=None, scores_pos_neg=None,
        name='', uselog=True, normed=True, default=-1, missing='?',
        linewidth=3, histtype='step', ncols=1, odds=False, **kwargs):
    if scores_pos_neg is None:
        pos_neg = pos_neg or pos_neg_from_arr(arr)
        scores_pos_neg = [x[name] for x in pos_neg]
    kwargs['bins'] = 30 if not 'bins' in kwargs else kwargs['bins']
    if odds:
        phist,nhist = [np.histogram(scores, range=(-1,1), density=True,
            **kwargs) for scores in scores_pos_neg]
        yvals = [np.nan_to_num(p/n) if n>-1 else -1 
                for p,n in zip(phist[0],nhist[0])]
        xvals = phist[1]
        plot(np.ravel(zip(xvals[:-1],xvals[1:])),
                np.ravel(zip(yvals,yvals)))
        return xvals,yvals
    else:
        (_,_,hp),(_,_,hn) = [hist(scores, log=uselog, histtype='step',
            linewidth=linewidth, normed=normed, **kwargs) 
            for scores in scores_pos_neg]
        return hp,hn

def presentation_mode(color='white', on=True):
    # customize individually with mpl.rcParams['text.color'] = '#000000'
    normmode = {
        'axes.facecolor': '#f8f8f8',
        'axes.labelcolor': '#222222',
        'xtick.color': '#222222',
        'ytick.color': '#222222',
        'figure.facecolor': '#dddddd',
        'axes.edgecolor': '#bcbcbc',
        'lines.linewidth': 2,
        #'axes.color_cycle': COLORSTRING, 
        #'text.color': '#222222'
        }
    whitemode = {
        'axes.facecolor': 'white',
        'axes.labelcolor': 'black',
        'xtick.color': '#000000',
        'ytick.color': '#000000',
        'figure.facecolor': 'white',
        'axes.edgecolor': 'black',
        'lines.linewidth': 3,
        #'axes.color_cycle': COLORSTRING, 
        #'text.color': '#222222'
        }
    blackmode = {
        'axes.facecolor': 'black',
        'axes.labelcolor': '#ffffff',
        'xtick.color': '#FFFFFF',
        'ytick.color': '#FFFFFF',
        'figure.facecolor': 'black',
        'axes.edgecolor': '#ffffff',
        'lines.linewidth': 2,
        #'axes.color_cycle': COLORSTRING_BLACK, 
        #'text.color': '#222222'
        }
    usemode = normmode if not on else (whitemode if color=='white' else
            blackmode)
    global COLORS
    COLORS = COLORS_BLACK if usemode == blackmode else COLORS_WHITE
    mpl.rcParams.update(usemode)
    
def ppis_scatter(ppis1, ppis2, useinds=range(3)):
    """
    useinds: set to [0,1,3,2] to take ppi.learning_examples output into (score,
    t/f) tuples; [0,1,3] to exclude the class.
    """
    pd1,pd2 = [pd.PairDict([[p[i] for i in useinds] for p in ppis]) 
            for ppis in ppis1,ppis2]
    nvals = len(useinds)-2
    pdcomb = pd.pd_union_disjoint_vals(pd1, pd2, adefaults=[0]*nvals,
            bdefaults=[0]*nvals)
    vals = zip(*ut.i1(pdcomb.d.items()))
    v1s,v2s = zip(*vals[:nvals]), zip(*vals[nvals:])
    v1s,v2s = [ut.i0(x) for x in v1s,v2s]
    return v1s,v2s

def scatter_union_labeled(avals, alabels, bvals, blabels):
    """
    vals are the columns of data to scatter (eg, el.mat[:,0]).  
    labels are el.prots.
    """
    dfs = [pa.DataFrame(data=vals,index=labels) for vals,labels in
            [(avals,alabels),(bvals,blabels)]]
    dfout = dfs[0].join(dfs[1], how='outer', rsuffix='_b')
    dfout = dfout.fillna(0)
    return dfout.values[:,0],dfout.values[:,1]

def eluts_scatter(elut1, elut2, ncols=None):
    vals1,vals2 = [],[]
    for icol in range(elut1.mat.shape[1])[:ncols]:
        new1, new2 = scatter_union_labeled(elut1.mat[:,icol], elut1.prots,
                elut2.mat[:,icol], elut2.prots)
        vals1 = np.concatenate((vals1, new1))
        vals2 = np.concatenate((vals2, new2))
    return vals1,vals2

def cluster_elut(mat):
    import hcluster
    ymat = hcluster.pdist(mat)
    zmat = hcluster.linkage(ymat)
    figure()
    order = hcluster.dendrogram(zmat)['leaves']
    clf() 
    imshow(mat[order,:])

def profiles_cxs(e, cxs, **kwargs):
    # blue/yellow/red map: 'jet'
    defaults = {'interpolation': 'nearest', 'cmap':'hot', 'vmin':1}
    kwargs = ut.dict_set_defaults(kwargs, defaults)
    arr = np.array(e.mat)
    dinds = ut.list_inv_to_dict(e.prots)
    useps = [p for c in cxs for p in c]
    useinds = [dinds[p] for p in useps if p in dinds]
    vals = np.clip(np.log2(arr[useinds,:]),0,100)
    imshow(vals, **kwargs)
    return vals

def scatter_blake(a, b, which='circles', classes=[0,1], colors=['k','r'],
        maxval=None, **kwargs):
    if maxval:
        a,b = filter_valpairs(a,b,maxval)
    defaults = {'s': 50, 'alpha':.2, 'lw':0}
    kwargs = ut.dict_set_defaults(kwargs, defaults)
    if type(a[0]) == list or type(a[0]) == tuple:
        # second value is presumed to be class--should be 0 or 1, which will be
        # mapped to the colormap cmap.
        # Also need to clean a and b to just be values rather than values and
        # classes.
        print 'using classes'
        assert ut.i1(a) == ut.i1(b), "Classes not the same between a and b"
        kwargs['c'] = [colors[0] if x==classes[0] else colors[1] for x in
                ut.i1(a)]
        a,b = ut.i0(a), ut.i0(b)
    else:
        c = 'k'
    if which=='pointcloud':
        scatter(a, b, s=50, alpha=0.08, lw=0)
        scatter(a, b, **kwargs)
    elif which=='points':
        scatter(a, b, **kwargs)
    elif which=='fadepoints':
        scatter(a, b, **kwargs)
    elif which=='circles':
        del kwargs['lw']
        scatter(a, b, facecolors='none', edgecolors=c, **kwargs)
    title('R-squared: %0.3f' % ut.r_squared(a,b))

def hexbin_blake(a, b, maxval=None, **kwargs):
    defaults = {'cmap': 'binary', 'bins':'log', 'gridsize':np.floor(sqrt(len(a))/2) }
    kwargs = ut.dict_set_defaults(kwargs, defaults)
    if maxval:
        a,b = filter_valpairs(a,b,maxval)
    hexbin(a, b, **kwargs)

def filter_valpairs(a,b,maxval):
    return zip(*[(x,y) for x,y in zip(a,b) if x<maxval and y<maxval])

def multi_scatter(comps,scatter_func=scatter_blake, preprocess=None,
        names=None, **kwargs):
    """
    Takes care of making subplots and labeling axes when comparing more than
    two sets of values.
    """
    total = len(comps)
    for i in range(total):
        for j in range(i+1,total):
            n = (total-1)*i+j
            print i,j,n
            subplot(total-1, total-1, n)
            ys,xs = comps[i],comps[j]
            # this syntax is mis-interpreted, and both new values go into xs
            #xs,ys = preprocess(xs,ys) if preprocess else xs,ys
            if preprocess:
                xs,ys = preprocess(xs,ys) 
            scatter_func(xs,ys, **kwargs)
            if names and j==i+1:
                ylabel(names[i])
                xlabel(names[j])

def hist_ndarray(arr, names=None, showindex=None, markerstyle='kp',
        do_hist=True, **kwargs):
    names = names or arr.dtype.names
    K = len(names)
    for i,name in enumerate(names):
        subplot(K,1,i+1)
        grid("off")
        if showindex is not None:
            plot(arr[name][showindex], 0, markerstyle, ms=15)
        if do_hist:
            hist(arr[name], **kwargs)
            legend([name], numpoints=1, loc=2)
            xlabel(name)
        #ylabel(name)

def hist_pairs_nonpairs(scores, pairs, negmult=1, do_plot=True, **kwargs):
    """
    scores: list of tuples (id1, id2, score)
    pairs: list of tuples (id1, id2)
    Make a histogram for scores of pairs against random sampling of non-pairs
    from the set of ids making up pairs.
    """
    assert len(pairs[0])==2, "Too many data points"
    nonpairs = pu.nonpairs_gen(pairs, len(pairs)*negmult)
    def scorelist_pairs(pairs, scores):
        pdscores = pd.PairDict([s[:3] for s in scores])
        for p in pairs:
            puse = pdscores.find(p)
            yield float(pdscores.d[puse][0]) if puse else 0
    pscores, nscores = [[x for x in scorelist_pairs(l, scores)] for l in pairs, nonpairs]
    if do_plot:
        hist(pscores, **kwargs)
        hist(nscores, **kwargs)
    return pscores, nscores

def fd_bincount(values):
    pass

def equal_freq_bins(values, nbins):
    pcts = arange(0, 100, 100/nbins)
    pctiles = [np.percentile(values, p) for p in pcts]
    return pctiles

def cdf(values, **kwargs):
    counts, edges = np.histogram(values, normed=True, **kwargs)
    counts = counts/sum(counts)
    cum_counts = np.cumsum(counts)
    plot(edges[1:], cum_counts)

def remove_nan(xs, ys):
    return zip(*[(x,y) for x,y in zip(xs,ys) if not(isnan(x)) and
        not(isnan(y))])

