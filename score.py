import sys
import numpy as np
import os
from scipy import sparse
import operator
import utils as ut
import elution as el

def score_examples(exstruct, score_mat, labels, name, default='?'):
    examples_out = []
    d = ut.list_inv_to_dict(labels)
    for e in exstruct.examples:
        p1 = e[0]
        p2 = e[1]
        if p1 in d and p2 in d:
            examples_out.append(e + [score_mat[d[p1],d[p2]]])
        else:
            examples_out.append(e + [default])
    exstruct.examples = examples_out
    exstruct.names.append(name)
    return exstruct

def scorekey_elution(score_key, elution, cutoff):
    if score_key == 'apex':
        score_mat = ApexScores(elution)
    elif score_key == 'poisson':
        score_mat = precalc_scores(elution.filename+'.corr_poisson', cutoff)
    elif score_key == 'wcc':
        score_mat = precalc_scores(elution.filename+'.T.wcc_width1', cutoff)
    else:
        assert False, "key not supported:" + score_key
    return score_mat
    
    
def score_examples_key(exstructs, score_key, elution, cutoff):
    score_mat = scorekey_elution(score_key, elution, cutoff)
    out = []
    for exstruct in exstructs:
        out.append(score_examples(exstruct, score_mat, elution.prots,
        score_key+'_'+ut.shortname(elution.filename)))
    return out 

def traver_corr(mat, repeat=200, norm='columns', verbose=True):
    # As described in supplementary information in paper.
    # Randomly draw from poisson(C=A+1/M) for each cell
    # where A = the observed count and M is the total fractions
    # normalize each column to sum to 1
    # then correlate, and average together for repeat tries.
    def poisson_corr(mat, iteration_display, norm):
        if verbose: print iteration_display
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


class ApexScores(object):

    def __init__(self, elution):
        mat = elution.mat
        self.apex_array = np.array(np.argmax(mat, axis=1))
        self.shape = (len(self.apex_array),len(self.apex_array))

    def __getitem__(self, index):
        return int(self.apex_array[index[0]] == self.apex_array[index[1]])

def precalc_scores(scoref, cutoff):
    save_sparse = ut.config()['save_sparse_corrs'] 
    sparsef = '%s.filt_%s.pyd' % (scoref, cutoff)
    if os.path.exists(sparsef): 
        return ut.loadpy(sparsef)
    else:
        mat = np.matrix(np.loadtxt(scoref))
        mat[mat < cutoff] = 0
        spmat = sparse.csr_matrix(mat)
        if save_sparse:
            print 'saving filtered ', sparsef
            ut.savepy(spmat, sparsef)
        return spmat


class CosineLazyScores(object):

    def __init__(self,elution):
        mat = elution.mat
        norms = np.apply_along_axis(np.linalg.norm, 1, mat)
        self.mat_rownormed = np.nan_to_num(mat / np.matrix(norms).T)
        assert type(self.mat_rownormed) == type(np.matrix(''))
        self.shape = (mat.shape[0],mat.shape[0])
        
    def __getitem__(self, index):
        # Dot product of normed rows
        return float(self.mat_rownormed[index[0],:] *
                    self.mat_rownormed[index[1],:].T)

    
def pairs_exceeding(elut, skey, thresh=0.5):
    if skey == 'apex':
        apex_obj = ApexScores(elut)
        apexes = apex_obj.apex_array
        pair_inds = [(proti,protj) for proti,peaki in enumerate(apexes) for
            protj,peakj in enumerate(apexes) if peaki==peakj and proti!=protj]
    else:
        # scorekey_elution now returns a csr sparse matrix
        # but conversion and where are quite fast anyway, so whatever
        score_mat = scorekey_elution(skey, elut, thresh).todense()
        rows, cols = np.where(score_mat > thresh)
        pair_inds =  ut.zip_exact(rows, cols)
    return pair_inds

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 3:
        sys.exit("usage: python corr.py filename method(poisson|dotproduct|corrcoef|cov) [argument]") 
    fname = sys.argv[1]
    method = sys.argv[2]
    methodarg = None if nargs < 4 else int(sys.argv[3])
    elut = el.load_elution(fname)
    if method == 'poisson':
        corr = traver_corr(elut.mat, repeat=methodarg) if methodarg else \
            traver_corr(elut.mat)
    elif method == 'dotproduct':
        corr = elut.mat * elut.mat.T
    elif method == 'corrcoef':
        corr = np.corrcoef(elut.mat)
    elif method == 'cov':
        corr = np.cov(elut.mat)
    fileout = fname+'.corr_'+method
    np.savetxt(fileout, corr, delimiter='\t')

