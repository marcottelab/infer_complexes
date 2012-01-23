import pylab
from pylab import *
import hcluster
import random

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
