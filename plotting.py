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
import hcluster
import pandas as pa


COLORS = ['#4571A8', 'black', '#A8423F', '#89A64E', '#6E548D', '#3D96AE',
           '#DB843D', '#91C4D5', '#CE8E8D', '#B6CA93', '#8EA5CB', 'yellow',
           'gray', 'blue']
          

def plot_result(result, ppis=None, **kwargs):
    ppis = ppis if ppis else result.ppis
    kwargs['label'] = kwargs.get('label','') +' '+ result.name
    pr_plot(ppis, result.ntest_pos, **kwargs)

def boot_resample(extr_exte):
    return [Struct(names=ex.names,examples=ut.sample_wr(ex.examples, len(ex.examples))) for ex in extr_exte]
    
def rolling_scores(tested, show=1000, window=50, **kwargs):
    #rolling = [len([t for t in tested[i:i+window] if t[3]==1])/window for i in
    #range(show-window)]
    padded = list(np.zeros((50,4)))+list(tested)
    rolling = [len([t for t in padded[i:i+window] if t[3]==1])/window for i in range(show)]
    #plot([0]*window + rolling, **kwargs)
    plot(rolling, **kwargs)
    plot([t[2] for t in tested[:show]], **kwargs)
    xlabel('starting index in scored examples')
    ylabel('fraction true in index:index+%s'%window)
    legend(['fraction true','score'])

def bar_plot(names, yvals, slanted_names=False, **kwargs):
    fig = figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(names)), yvals, align='center', **kwargs)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    if slanted_names: fig.autofmt_xdate()
    show()
    return ax
    
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
    kwargs['label'] = kwargs.get('label','') + ' %.3f' % cv.auroc(xs,ys)
    plot(xs, ys, **kwargs)
    plot([0,xs[-1]], [0,ys[-1]], 'k--')
    
def pr_plot(cv_pairs, total_trues, rescale=None, prec_test=None, **kwargs):
    """
    rescale: adjust precision values assuming rescale times as many negatives
    total_trues:
    - None for just displaying recall count instead of fraction
    - 'auto' to calculate from the supplied tested cv_pairs
    - integer to use that supplied integer as total trues
    """
    if total_trues == 'auto':
        total_trues = len([t for t in cv_pairs if t[3]==1])
    recall,precision = cv.pr(cv_pairs) 
    if rescale:
        precision = [ p / (p + (1-p) * rescale) for p in precision]
    if prec_test:
        kwargs['label'] = kwargs.get('label','') + (' Re:%0.2f' %
        cv.calc_recall(precision,prec_test, total_trues)) + (' @ Pr:%0.2f'
            % precision_test)
    if total_trues:
        recall = [r/total_trues for r in recall]
    plot(recall, precision, **kwargs)
    xlabel('Recalled Correctly')
    ylabel('Precision: TP/(TP+FP)')
    ylim(-0.01,1.01)
    xlim(xmin=-0.002)
    legend()

def imshow2(*args):
    imshow(*args, interpolation='nearest', aspect='auto',
           cmap='bone', vmin=0)
           # Need vmin = 0 so the lowest values aren't represented 
           # as say black by the colorbar if they're not 0.
           # Colormaps: bone, gray

def examples_dist_arr(arr, score_indices, uselog=True, normed=True, default=-1,
                   missing='?', linewidth=3, histtype='step', ncols=1,
                   **kwargs):
    nplots = len(score_indices)+1
    pos,neg = [arr[[i for i in range(len(arr)) if arr[i][2]==t]] for t in 1,0]
    for i, (ind, name) in enumerate([(ind,arr.dtype.names[ind]) for ind in score_indices]):
        subplot(int(nplots/ncols)+1, ncols, i+2)
        #kwargs['range'] = kwargs['range'] if 'range' in kwargs else \ [func([func(data) for data in [pos,neg]]) for func in [min,max]]
        kwargs['bins'] = 30 if not 'bins' in kwargs else kwargs['bins']
        (_,_,hp),(_,_,hn) = [hist(data, log=uselog, histtype='step', linewidth=linewidth, normed=normed, **kwargs) for data in pos[name],neg[name]]
        title(name)
    subplot(int(nplots/ncols)+2,ncols,1)
    legend([hp[0],hn[0]],['Pos','Neg'])
        
def presentation_mode(on=True):
    # customize individually with mpl.rcParams['text.color'] = '#000000'
    normmode = {
        'axes.facecolor': '#f8f8f8',
        'axes.labelcolor': '#222222',
        'figure.facecolor': '#dddddd',
        'axes.edgecolor': '#bcbcbc',
        'lines.linewidth': 2,
        #'text.color': '#222222'
        }
    presmode = {
        'axes.facecolor': 'white',
        'axes.labelcolor': '#000000',
        'figure.facecolor': 'white',
        'axes.edgecolor': '#000000',
        'lines.linewidth': 3,
        #'text.color': '#222222'
        }
    usemode = presmode if on else normmode
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
    return v1s,v2s

def eluts_scatter(avals, alabels, bvals, blabels):
    """
    vals are the columns of data to scatter (eg, el.mat[:,0]).  
    labels are el.prots.
    """
    dfs = [pa.DataFrame(data=vals,index=labels) for vals,labels in
            [(avals,alabels),(bvals,blabels)]]
    dfout = dfs[0].join(dfs[1], how='outer', rsuffix='_b')
    dfout = dfout.fillna(0)
    return dfout.values[:,0],dfout.values[:,1]

def cluster_elut(mat):
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

def scatter_blake(a, b, which='fadepoints', classes=[0,1], colors=['k','r'],
        **kwargs):
    defaults = {'s': 4, 'lw':0}
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
        scatter(a, b, s=50, alpha=0.2, lw=0)
    title('R-squared: %0.3f' % ut.r_squared(a,b))

def hexbin_blake(a, b, **kwargs):
    defaults = {'cmap': 'binary', 'bins':'log'}
    kwargs = ut.dict_set_defaults(kwargs, defaults)
    hexbin(a, b, gridsize=np.floor(sqrt(len(a))/2), **kwargs)

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


