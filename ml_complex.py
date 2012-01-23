from __future__ import division
import numpy as np
import operator

def coelution(elut, norm=1):
    # elut.mat: float matrix of spectral counts
    # elut.genes: gene labels for the rows of the array
    if norm:
        # stupid simple: pearson correlation matrix
        corr = np.corrcoef(elut.mat) # between -1 and 1
    else:
        corr = np.cov(elut.mat)
    return corr

def traver_corr(mat, repeat=1000):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each row to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_corr(mat):
        M = mat.shape[1]
        C = mat + 1/M
        poisson_mat = np.matrix(np.zeros(C.shape))
        for i in range(C.shape[0]):
            for j in range(M):
                poisson_mat[i,j] = np.random.poisson(C[i,j])
        normp_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
        corr = np.nan_to_num(np.corrcoef(normp_mat))
        return corr
    avg_result = (reduce(operator.add, (poisson_corr(mat) for i in range(repeat)))
               / repeat)
    return avg_result
