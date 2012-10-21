from __future__ import division
import numpy as np
import pylab
from pylab import *
import random
import itertools
import cv
from Struct import Struct
import utils as ut
import ppi
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
    
