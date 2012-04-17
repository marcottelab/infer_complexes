from __future__ import division
import numpy as np

def cv_pairs(scores, true_pairs, genes, sample_frac=.1):
    # scores: input 2d array of scores for each index-index pair
    # true_pairs: dict: {gene1: set(gene2,gene3,gene4), gene2:
    #   set(gene1,gene5)}
    ranked = rank_scores(scores, sample_frac)
    tested = []
    for i,j,score in ranked:
        if i!=j: # remove same-same
            hit = int(genes[j] in true_pairs.get(genes[i],[]))
            tested.append( (genes[i], genes[j], score, hit) )
    return tested

def rank_scores(scores, sample_frac=.1):
    print "cv.rank_scores(): Should Redo if only top N are desired."
    import random
    ranked_scores = []
    for i in range(scores.shape[0]):
        for j in range(scores.shape[1]):
            ranked_scores.append((i,j,scores[i,j]))
    ranked_scores = random.sample(ranked_scores, int(scores.size * sample_frac))
    ranked_scores.sort(key = lambda x: x[2], reverse = True)
    return ranked_scores

def roc(tested_pairs):
    # tested pairs: [ (row, col, score, hit(0/1)), ...]
    x = 0
    y = 0
    xs = [x]
    ys = [y]
    score_prev = 0
    for _,_,score,hit in tested_pairs:
        if hit==1:
            y += 1
        else:
            x += 1
        # We have this conditional so that segments with equal probability
        # are drawn as a line with a + at the end.
        if score != score_prev:
            xs.append(x)
            ys.append(y)
            score_prev = score
    # don't forget to add the last point
    xs.append(x)
    ys.append(y)
    return xs,ys

def pr(tested_pairs):
    # precision-recall
    # tested pairs: [ (row, col, score, hit(0/1)), ...]
    hits = np.array([tp[3] for tp in tested_pairs])
    hit_inds = np.where(hits==1)[0]+1
    precision = [(i+1)/(tp_plus_fp) for i,tp_plus_fp in
                 enumerate(hit_inds)]
    return range(1,len(precision)+1), precision

def examples_to_cvpairs(exlist, scoreindex=None):
    """
    Translates and reorders an example list into the format used in cv
    columns to give (id1, id2, score, hit(1/0)).  And sorts.
    """
    def hit_translate(tf):
        return 1 if tf == 'true' else 0
    scoreindex = scoreindex if scoreindex else len(exlist.examples[0])-1
    print "Sample score: %s", exlist.examples[0][scoreindex]
    reworked_examples = [(e[0],e[1],e[scoreindex],hit_translate(e[2])) for e in
        exlist.examples]
    reworked_examples.sort(key=lambda x:x[2],reverse=True)
    return reworked_examples

def calc_recall(precisions, gt_value):
    passing_inds = np.where( np.array(precisions) >= gt_value )[0]
    return np.max(passing_inds) if len(passing_inds)>0 else 0
    
def auroc(xs,ys):
    auroc = 0
    xprev = 0
    yprev = 0
    for x,y in zip(xs,ys):
        auroc += .5*(y+yprev)*(x-xprev)
        xprev = x
        yprev = y
    auroc = auroc / (x*y)
    return auroc

        
