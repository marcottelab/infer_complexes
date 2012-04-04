from __future__ import division
import pylab
from pylab import *
import random
import hcluster
import cv
COLORS = ['#4571A8', '#A8423F', '#89A64E', '#6E548D', '#3D96AE', '#DB843D',
           '#91C4D5', '#CE8E8D', '#B6CA93', '#8EA5CB', 'yellow', 'gray',
           'blue', 'black']

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
    
def pr_plot(cv_pairs, precision_test, **kwargs):
    precision_test = precision_test if precision_test else 0.95
    recall,precision = cv.pr(cv_pairs) 
    kwargs['label'] = kwargs.get('label','') + (' Re:%i' %
        cv.calc_recall(precision,precision_test))
    plot(recall, precision, **kwargs)

def roc_plot_examples(exlist, scoreindex, **kwargs):
    roc_plot(cv.examples_to_cvpairs(exlist, scoreindex), **kwargs)
    
def pr_plot_examples(exlist, scoreindex, precision_check=.5, **kwargs):
    pr_plot(cv.examples_to_cvpairs(exlist, scoreindex), precision_check,
        **kwargs)

def imshow2(*args):
    imshow(*args, interpolation='nearest', aspect='auto',
           cmap='bone', vmin=0)
           # Need vmin = 0 so the lowest values aren't represented 
           # as say black by the colorbar if they're not 0.
           # Colormaps: bone, gray

def examples_dist(exlist,score_indices):
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
    
