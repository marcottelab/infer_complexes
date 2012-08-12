import sys
import itertools
import numpy as np
import os
from scipy import sparse
from collections import defaultdict
import operator
import utils as ut
import elution as el
import orth


def score_array_multi(arr, sp_base, elut_fs, scores, cutoff, verbose=False):
    eluts = [(el.load_elution(f),f) for f in elut_fs]
    print "Fix 2, try2."
    current_sp = ''
    for e,f in eluts:
        new_sp = os.path.basename(f)[:2]
        if new_sp != current_sp:
            print "Starting first %s file: %s" % (new_sp, os.path.basename(f))
            current_sp = new_sp
        sp_target = ut.shortname(f)[:2]
        baseid2inds = orth_indices(sp_base, sp_target, e.prots)
        for score in scores:
            if verbose: print score, f
            score_array(arr, e, f, score, cutoff, baseid2inds)

def orth_indices(sp_base, sp_target, prot_list):
    """
    Using appropriate orthology, take a list of target species gene ids
    (corresponding to rows in the target species score matrix), and
    return a dict mapping base species gene ids to (sets of) indices in that
    list and therefore to (sets of) row/column indices in the square
    interaction score matrix. 
    """
    targ2inds = dict([(k,set([v]))
                      for k,v in ut.list_inv_to_dict(prot_list).items()])
    if sp_base == sp_target:
        return targ2inds
    else:
        base2targ = orth.odict(sp_base, sp_target)
        base2inds = ut.compose_dict_sets(base2targ, targ2inds)
        base2inds = dict([(k,v) for k,v in base2inds.items() if len(v)>0])
        return base2inds

def score_array(arr, elut, fname, score, cutoff, id2inds):
    """
    Use the target species score matrix to get interaction pair in the base
    species array.  Don't score and just leave as default (0 now) cases where
    either: 1) One of the pair is not in this score matrix, or 2) The two base
    ids in the pair map to identical targets, since in that case we also can
    get no information from this data (see notes 2012.08.12).
    """
    if score == 'apex':
        score_mat = ApexScores(elut)
    else:
        fscore = fname + (
                  '.corr_poisson' if score=='poisson' else
                  '.T.wcc_width1' if score=='wcc' else
                  0 ) # no score: exception since string and int don't add
        score_mat = precalc_scores(fscore)
    score_name = name_score(fname,score)
    for i,row in enumerate(arr):
        id1,id2 = row['id1'],row['id2']
        #if id1 in id2inds and id2 in id2inds and id2inds[id1]!=id2inds[id2]: 
        if id1 in id2inds and id2 in id2inds: 
            # Could also check for i!=j but would have no effect here since
            # these mappings come from disjoint orthogroups.
            #row[score_name] = max([score_mat[i,j] for i in id2inds[id1] for j in id2inds[id2]])
            scores = [score_mat[i,j] for i in id2inds[id1] for j in id2inds[id2] if i!=j]
            if len(scores)>0: 
                row[score_name] = max(scores)
        
def name_score(fname, score):
    return ut.shortname(fname) + '_' + score 

def scorekey_elution(score_key, elution):
    if score_key == 'apex':
        score_mat = ApexScores(elution)
    elif score_key == 'poisson':
        score_mat = precalc_scores(elution.filename+'.corr_poisson')
    elif score_key == 'wcc':
        score_mat = precalc_scores(elution.filename+'.T.wcc_width1')
    else:
        assert False, "key not supported:" + score_key
    return score_mat
    
    
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

def precalc_scores(scoref, dtype='f2'):
    # NOTE to change dtype you must change it in loadtxt below!!
    save_compact = ut.config()['save_compact_corrs'] 
    compactf = '%s.%s.pyd' % (scoref, dtype)
    if os.path.exists(compactf): 
        return ut.loadpy(compactf)
    else:
        ascores = np.loadtxt(scoref, dtype='f2')
        if save_compact:
            print 'saving compact', compactf
            ut.savepy(ascores, compactf)
        return ascores


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
    else: # loading precomputed indices is so far massively slower than this
        score_mat = scorekey_elution(skey, elut)
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

