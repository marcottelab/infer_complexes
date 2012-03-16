import pylab
from pylab import *
import random
import hcluster
import cv

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

def roc_plot(tested_pairs, **kwargs):
    xs,ys = cv.roc(tested_pairs)
    plot(xs, ys, label=' %.3f' % cv.auroc(xs,ys), **kwargs)
    plot([0,xs[-1]], [0,ys[-1]], '--')

def imshow2(*args):
    imshow(*args, interpolation='nearest', aspect='auto',
           cmap='bone', vmin=0)
           # Need vmin = 0 so the lowest values aren't represented 
           # as say black by the colorbar if they're not 0.
           # Colormaps: bone, gray
