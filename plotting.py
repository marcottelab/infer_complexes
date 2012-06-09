from __future__ import division
import pylab
from pylab import *
import random
import itertools
import hcluster
import cv
from Struct import Struct
import utils as ut
import myml
import ppi
COLORS = ['#4571A8', '#A8423F', '#89A64E', '#6E548D', '#3D96AE', '#DB843D',
           '#91C4D5', '#CE8E8D', '#B6CA93', '#8EA5CB', 'yellow', 'gray',
           'blue', 'black']

def pr_ppi(train_test, ntest_pos, prec_check=0, label_stats=True,
           cutoff_score=None, **kwargs):
    tested = myml.fit_and_test(train_test)
    kwargs['label'] = kwargs.get('label','')+" %s Total Test Pos; " % \
           ntest_pos + ppi.exstats(extr_exte)
    pr_plot(tested, prec_check, ntest_pos, label_prec=False, **kwargs)

def boot_resample(extr_exte):
    return [Struct(names=ex.names,examples=ut.sample_wr(ex.examples, len(ex.examples))) for ex in extr_exte]
    
def score_threshold(tested, show=1000, window=50):
    rolling = [len([t for t in tested[i:i+window] if t[3]==1])/window for i in range(show-window)]
    plot([0]*window + rolling)
    plot([t[2] for t in tested[:show]])
    xlabel('starting index in scored examples')
    ylabel('fraction true in index:index+%s'%window)
    
def cluster(corr):
    # corr: a matrix of similarity scores, such as a covariance matrix
    ymat = hcluster.pdist(corr)
    zmat = hcluster.linkage(ymat)
    figure()
    order = hcluster.dendrogram(zmat)['leaves']
    figure()
    imshow(corr[order,:][:,order])
    # check for failure signs
    for i in random.sample(range(len(order)),10):
        if order[i] - order[i-1] == 1:
            print 'HEY!! probable clustering failure.'
            break
    return order

def roc_plot(cvpairs, **kwargs):
    xs,ys = cv.roc(cvpairs) 
    kwargs['label'] = kwargs.get('label','') + ' %.3f' % cv.auroc(xs,ys)
    plot(xs, ys, **kwargs)
    plot([0,xs[-1]], [0,ys[-1]], 'k--')
    
def pr_plot(cv_pairs, precision_test, total_trues, rescale=None, label_prec=True, **kwargs):
    """
    rescale: adjust precision values assuming rescale times as many negatives
    total_trues:
    - None for just displaying recall count instead of fraction
    - 'auto' to calculate from the supplied tested cv_pairs
    - integer to use that supplied integer as total trues
    """
    precision_test = precision_test if precision_test else 0.0
    if total_trues == 'auto':
        total_trues = len([t for t in cv_pairs if t[3]==1])
    recall,precision = cv.pr(cv_pairs) 
    if rescale is not None:
        precision = [ p / (p + (1-p) * rescale) for p in precision]
    if label_prec:
        kwargs['label'] = kwargs.get('label','') + (' Re:%0.2f' %
        cv.calc_recall(precision,precision_test, total_trues)) + (' @ Pr:%0.2f'
            % precision_test)
    if total_trues is not None:
        recall = [r/total_trues for r in recall]
    plot(recall, precision, **kwargs)
    xlabel('Recalled Correctly')
    ylabel('Precision: TP/(TP+FP)')
    ylim(0,1.02)

def pr_plot_weka(fname, prec_test=0.0, total_trues=None, rescale=None):
    pr_plot(cv.load_weka_filtered_tpairs(fname), prec_test, total_trues,
        rescale=rescale, label=ut.shortname(fname)) 
    
def roc_plot_examples(exlist, scoreindex, **kwargs):
    roc_plot(cv.examples_to_cvpairs(exlist, scoreindex), **kwargs)
    
def pr_plot_examples(exlist, scoreindex, precision_check=.5, plot_random=True, **kwargs):
    cvs = cv.examples_to_cvpairs(exlist, scoreindex)
    pr_plot(cvs, precision_check, **kwargs)
    if plot_random:
        random.shuffle(cvs)
        pr_plot(cvs, None, color='k', linestyle='--', label='shuffled')

def imshow2(*args):
    imshow(*args, interpolation='nearest', aspect='auto',
           cmap='bone', vmin=0)
           # Need vmin = 0 so the lowest values aren't represented 
           # as say black by the colorbar if they're not 0.
           # Colormaps: bone, gray

def examples_dist(exlist, score_indices, uselog=True, normed=True, default=-1,
                   missing='?', linewidth=3, histtype='step', **kwargs):
    hs = []
    nplots = len(score_indices)
    for i, (ind, name) in enumerate([(ind,exlist.names[ind]) for ind in score_indices]):
        subplot(nplots,1,i+1)
        pos,neg = [[float(e[ind]) if e[ind]!=missing else default for e in
                exlist.examples if e[2]==truth] for truth in ['true','false']]
        kwargs['range'] = kwargs['range'] if 'range' in kwargs else \
                   [func([func(data) for data in [pos,neg]]) for func in
                   [min,max]]
        kwargs['bins'] = 30 if not 'bins' in kwargs else kwargs['bins']
        (_,_,hp),(_,_,hn) = [hist(data, log=uselog, histtype='step',
                   linewidth=linewidth, normed=normed, **kwargs) for data in
                   zip([pos,neg])]
        legend([hp[0],hn[0]],['Pos','Neg'],loc=3)
        title(name)
        
def examples_dist_old(exlist,score_indices):
    hs = []
    nplots = len(score_indices)
    for i, (ind, name) in enumerate([(ind,exlist.names[ind-3]) for ind in score_indices]):
        subplot(nplots,2,2*i+1)
        pos = [float(e[ind]) for e in exlist.examples if e[2]=='true' and
                e[ind]!='?']
        neg = [float(e[ind]) for e in exlist.examples if e[2]=='false' and
    e[ind]!='?']
        _,_,hp = hist(pos, linewidth=0, bins=30, log=True, alpha=0.5)
        _,_,hn = hist(neg, linewidth=0, bins=30, log=True, alpha=0.5)
        legend([hp[0],hn[0]],['Pos','Neg'])
        pct_excluded = 100 * (1 - ((len(pos) + len(neg)) / len(exlist.examples)))
        title(name+' log, not stacked, %ipct exs excluded' % pct_excluded)
        subplot(nplots,2,2*i)
        _,_,h = hist((pos,neg), linewidth=0, histtype='barstacked', bins=30)
        legend([h[0][0],h[1][0]],['Pos','Neg'])
        title(name+' linear, stacked')

def examples_scatter_posneg(exlist, baseindex, otherindices):
    """
    Scatter the scores in the examplelist.  Scatter all scores against
    baseindex. First score starts at index 3.
    """
    exs = exlist.examples
    first_index = 3
    ncols = len(exs[0])
    for i,index in enumerate([i for i in otherindices if i!=baseindex]):
        subplot(len(otherindices),1,i+1)
        (pos, negs) = filter_posnegs(exs, baseindex, index)
        r2 = np.corrcoef(*zip(*(pos+negs)))[0][1]
        # s controls marker size, lw0 gets rid of any surrounding line
        scatter(*zip(*negs),color='r',s=2,lw=0)
        scatter(*zip(*pos),color='b',s=2,lw=0)
        xlabel(exlist.names[baseindex-first_index])
        ylabel(exlist.names[index-first_index])
        title('R2 = %0.2f' % r2)
        legend(['R2 = %0.2f' % r2])
            
def examples_scatter_scores(exlist, baseindex=3):
    """
    Scatter the scores in the examplelist.  Scatter all scores against
    baseindex. First score starts at index 3.
    """
    exs = exlist.examples
    first_index = 3
    ncols = len(exs[0])
    for i,index in enumerate([i for i in range(first_index,ncols) if i!=baseindex]):
        filtered = filter_unknowns(exs, baseindex, index)
        r2 = np.corrcoef(*zip(*filtered))[0][1]
        scatter(*zip(*filtered), color=COLORS[i],
            label=exlist.names[index-first_index]+'; R2=%0.2f' % r2,alpha=0.5)
    title('Scatter against ' + exlist.names[baseindex-first_index])
    xlabel(exlist.names[baseindex-first_index])
    legend()
            
def filter_posnegs(exs, i1, i2, unknown='?'):
    negs = [(e[i1],e[i2]) for e in exs if e[i1]!=unknown and e[i2]!=unknown and e[2]=='false']
    pos =  [(e[i1],e[i2]) for e in exs if e[i1]!=unknown and e[i2]!=unknown and e[2]=='true']
    return pos, negs
    
def filter_unknowns(exs, i1, i2, unknown='?'):
    return [(e[i1],e[i2]) for e in exs if e[i1]!=unknown and e[i2]!=unknown] 

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
    
