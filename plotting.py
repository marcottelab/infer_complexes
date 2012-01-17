import pylab
from pylab import *
import hcluster

def cluster(corr):
    # corr: a matrix of similarity scores, such as a covariance matrix
    ymat = hcluster.pdist(corr)
    zmat = hcluster.linkage(ymat)
    figure()
    order = hcluster.dendrogram(zmat)['leaves']
    figure()
    imshow(corr[order,:][:,order])
