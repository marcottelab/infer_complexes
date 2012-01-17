import numpy as np

def coelution(elut, norm=1):
    # elut.arr: float array of spectral counts
    # elut.genes: gene labels for the rows of the array
    if norm:
        # stupid simple: pearson correlation matrix
        corr = np.corrcoef(elut.arr) # between -1 and 1
    else:
        corr = np.cov(elut.arr)
    return corr
