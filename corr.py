import numpy as np

def coelution(elut, norm=1):
    # elut.mat: float matrix of spectral counts
    # elut.genes: gene labels for the rows of the array
    if norm:
        # stupid simple: pearson correlation matrix
        corr = np.corrcoef(elut.mat) # between -1 and 1
    else:
        corr = np.cov(elut.mat)
    return corr

def traver_corr(mat, repeat=200, norm='columns'):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each column to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_corr(mat, iteration_display, norm=True):
        print iteration_display
        M = mat.shape[1]
        C = mat + 1/M
        poisson_mat = np.matrix(np.zeros(C.shape))
        for i in range(C.shape[0]):
            for j in range(M):
                poisson_mat[i,j] = np.random.poisson(C[i,j])
        if norm=='columns': 
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 0))
        elif norm=='rows': # seems to make no performance difference 1/25
            poisson_mat = np.nan_to_num(poisson_mat / np.sum(poisson_mat, 1))
        corr = np.nan_to_num(np.corrcoef(poisson_mat))
        return corr
    avg_result = (reduce(operator.add, (poisson_corr(mat, i, norm=norm) for i in
                                        range(repeat))) / repeat)
    return avg_result

