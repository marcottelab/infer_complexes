import sys
import itertools
import numpy as np
import os
from scipy import sparse
from collections import defaultdict
import operator
import utils as ut
import elution as el


def scores_array(arr, elut_fs, scores, cutoff):
    eluts = [(el.load_elution(f),f) for f in elut_fs]
    for e,f in eluts:
        e.idict = ut.list_inv_to_dict(e.prots)
        for score in scores:
            print score, f
            score_array_precomp(arr, e, f, score, cutoff)

def score_array_precomp(arr, elut, fname, score, cutoff):
    if score == 'apex':
        score_mat = ApexScores(elut)
    else:
        fscore = fname + (
                  '.corr_poisson' if score=='poisson' else
                  '.T.wcc_width1' if score=='wcc' else
                  0 ) # no score: exception since string and int don't add
        score_mat = precalc_scores(fscore, cutoff)
    idict = elut.idict
    name = name_score(fname,score)
    for i,row in enumerate(arr):
        p1,p2 = row['id1'],row['id2']
        if p1 in idict and p2 in idict:
            row[name] = score_mat[idict[p1],idict[p2]]
        
def name_score(fname, score):
    return ut.shortname(fname) + '_' + score 

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
                name_score(score_key, elution.filename)))
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
        self.apex_array = np.argmax(np.array(elution.mat), axis=1)
        self.shape = (len(self.apex_array),len(self.apex_array))

    def __getitem__(self, index):
        return int(self.apex_array[index[0]] == self.apex_array[index[1]])

def precalc_scores(scoref, cutoff):
    save_sparse = ut.config()['save_sparse_corrs'] 
    sparsef = '%s.filt_%s.pyd' % (scoref, cutoff)
    print cutoff, sparsef
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

def matching_pairs(values):
    """
    Return all pairs of indices in the given list whose values match
    """
    d = defaultdict(list)
    for ind,val in enumerate(values):
        d[val].append(ind)
    return [(i,j) for value in d for i,j in itertools.combinations(d[value],2)]
    
def pairs_exceeding(elut, skey, thresh):
    if skey == 'apex':
        apexes = ApexScores(elut).apex_array
        pair_inds = matching_pairs(apexes)
    else:
        # scorekey_elution now returns a csr sparse matrix
        rows, cols = scorekey_elution(skey, elut, thresh).nonzero()
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

