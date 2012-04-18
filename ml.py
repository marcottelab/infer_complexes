import random
import os
import utils as ut
from Struct import Struct
import complex as co
import elution as el
import fnet

def shuffled_base(positives, negatives):
    """
    positives, negatives: like [('id1','id2','true'/'false'), ...]
    """
    examples = positives + negatives
    random.shuffle(examples)
    exstruct = Struct(examples=examples, names=['id1','id2','hit'])
    return exstruct

def _dep_examples_from_scores(positives, negatives, mats_labels_names):
    """
    Generate a list of [id1, id2, hit, score1, score2, ...]
    positives: list of true interaction pairs
    mats_labels: list of tuples of [(score_matrix, row_labels, score_name), ..]
    Kneg: number of negative examples to generate from items in the pairs
    example:
    ex1 =
    ml.interaction_examples(random.sample(corumpairsen,1000),[(ellong.corr,
    ellong.prots,'longest'), (elall.corr, elall.prots, 'all')], 1000)
    """
    examples = []
    negatives = negatives if negatives else negative_pairs(positives,len(positives))
    def unzip(pairlist):
        return [p[0] for p in pairlist], [p[1] for p in pairlist]
    for (pairs,iscomplex) in ((positives, ['true']*len(positives)), (negatives,
        ['false']*len(negatives))):
        scores = [scores_interaction_matrix(pairs, mat, labels) for
                    (mat, labels, name) in mats_labels_names]
        p0,p1 = unzip(pairs)
        examples += [list(a) for a in zip(p0, p1, iscomplex, *scores)]
    # shuffling is necessary, otherwise ordering persists when scores are the
    # same.  see data 3/29.
    random.shuffle(examples) 
    example_list = Struct(examples=examples, names=['id1','id2','hit']+[mln[2]
                for mln in mats_labels_names])
    return example_list

def _dep_scores_interaction_matrix(pairs, mat, labels, default='?'):
    # pairs: interacting partners: [(a,b),(a,d),(b,c),...]
    # mat: N by N matrix of interaction scores
    # labels: length N list of labels for rows/columns (same)
    assert mat.shape[0] == mat.shape[1], "score matrix must be symmetric"
    labeldict = ut.list_inv_to_dict(labels)
    scores = []
    for a,b in pairs:
        if a in labels and b in labels:
            scores.append(mat[labeldict[a],labeldict[b]])
        else:
            scores.append(default)
    return scores

def examples_combine_scores(exlist, index_start, index_end, reduce_func,
        retain_scores=False, default=0, unknown='?'):
    """
    Combines scores at the specified index columns using reduce_func.
    Put the new score as a new column on the end or replaces with new.
    """
    exs = exlist.examples
    def replace_unknowns(lst):
        newlist = []
        for l in lst:
            newval = l if l!=unknown else default
            newlist.append(newval)
        return newlist
    newexs = [e if retain_scores else e[:index_start] + [reduce(reduce_func,
        replace_unknowns(e[index_start:index_end]))] for e in exs]
    names = (exlist.names if retain_scores else exlist.names[:index_start]) + \
            [exlist.names[index_start] + reduce_func.__name__]
    newex_struct = Struct(examples=newexs, names=names)
    return newex_struct

def weka_export(exstruct, filename, startindex=3, howmany=None):
    # skip the first 2 items specifiying the interaction
    # Note that I"m taking iscomplex(true/false) and putting it at the end
    # per weka convention
    f = open(filename,'w')
    f.write('@RELATION complexes\n\n')
    for name in exstruct.names[startindex:]:
        f.write('@ATTRIBUTE '+name+'\treal\n')
    f.write('@ATTRIBUTE iscomplex\t{true,false}\n')
    f.write('\n@DATA\n\n')
    for ex in exstruct.examples[:howmany]:
        f.write(', '.join([str(val) for val in ex[3:]]+[ex[2]])+'\n')
    f.write('%\n%\n%\n')
    f.close()


def _dep_negative_pairs(truepairs, K):
    # TODO: speed up duplicate checking
    items = list(reduce(set.union,[set(p) for p in truepairs]))
    negatives = []
    dtrues = pairs_to_dict(truepairs)
    dnegs = {}
    def pair_in(pair, dict):
        return ( pair[0] in dict.get(pair[1],set([])) or pair[1] in
                 dict.get(pair[0],set([])) )
    while len(negatives) < K:
        a = random.choice(items) 
        b = random.choice(items)
        newpair = (a,b)
        if a!=b and not ( pair_in(newpair, dtrues) or
                          pair_in(newpair, dnegs)):
            negatives.append(newpair)
            dnegs.setdefault(a,set([])).add(b)
    return negatives

def pairs_to_dict(pairs):
    d = {}
    for i,j in pairs:
        d.setdefault(i,set([])).add(j)
    return d

def combined_examples(pos_pairs, negatives, elutions, score_func, combine_func,
    retain_scores=False):
    col1 = 3
    combine_col = 3 + len(elutions)
    examples = ml.examples_from_scores(pos_pairs, negatives,
        [(score_func(e), e.prots, e.name) for e in elutions])
    return ml.examples_combine_scores(examples, col1, combine_col,
        combine_func, retain_scores=retain_scores)

def full_examples(poskey, npos, negskey, nnegs, elut_fs, scores, entrez_dict,
                  out_base='weka/'):
    pos = co.pairs_key(poskey)[:npos]
    if nnegs==None and npos==None: nnegs=len(pos)
    negs = co.pairs_key(negskey)[:nnegs]
    ex_struct = shuffled_base(pos,negs) 
    el.score_multi_elfs(ex_struct, elut_fs, scores)
    fnet.score_examples(ex_struct, genedict=entrez_dict)
    out_fname = os.path.join(out_base,
                  poskey+str(npos)+'_'+str(nnegs)+'negs.arff')
    if os.path.exists(out_fname):
        out_fname = ut.pre_ext(out_fname,str(random.randint(0,100)))
    weka_export(ex_struct, out_fname)
    print 'ready:', out_fname
