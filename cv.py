from __future__ import division
import numpy as np
import sys
import utils as ut
import pairdict as pd

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

def roc(tested_pairs, total_pos_neg=None):
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
    if total_pos_neg:
        xs.append(total_pos_neg[0])
        ys.append(total_pos_neg[1])
    return xs,ys

def pr(tested_pairs):
    # precision-recall
    # tested pairs: [ (row, col, score, hit(0/1)), ...]
    hits = np.array([tp[3] for tp in tested_pairs])
    hit_inds = np.where(hits==1)[0]+1
    precision = [(i+1)/(tp_plus_fp) for i,tp_plus_fp in
                 enumerate(hit_inds)]
    return range(1,len(precision)+1), precision

def aupr(tested, ntest_pos):
    """
    Excludes calculating the tail. Should neglibly affect calculations,
    especially in terms of comparisons, although worth noting it's wrong.
    """
    recall, precision = pr(tested)
    undercurve = sum(precision) / ntest_pos
    return undercurve

def examples_to_cvpairs(exlist, scoreindex=None):
    """
    Translates and reorders an example list into the format used in cv
    columns to give (id1, id2, score, hit(1/0)).  And sorts.
    """
    scoreindex = scoreindex if scoreindex else len(exlist.examples[0])-1
    print "Sample score: %s", exlist.examples[0][scoreindex]
    reworked_examples = [(e[0],e[1],e[scoreindex],true_to_1(e[2])) for e in
        exlist.examples]
    reworked_examples.sort(key=lambda x:x[2],reverse=True)
    return reworked_examples

def gold_label_ppis(ppis, gold_ppis):
    pdgold = pd.PairDict(gold_ppis)
    return [(p[0],p[1],p[2],1 if pdgold.contains((p[0],p[1])) else 0) for p in
            ppis]

def true_to_1(tf):
    return 1 if tf == 'true' else 0

def calc_recall(precisions, gt_value, recall_rate_of_total=None):
    passing_inds = np.where( np.array(precisions) >= gt_value )[0]
    recalled = np.max(passing_inds) if len(passing_inds)>0 else 0
    if recall_rate_of_total is not None:
        recalled /= recall_rate_of_total
    return recalled
    
def auroc(xs,ys):
    if len(xs)==0: return 0
    auroc = 0
    xprev = 0
    yprev = 0
    for x,y in zip(xs,ys):
        auroc += .5*(y+yprev)*(x-xprev)
        xprev = x
        yprev = y
    auroc = auroc / (x*y)
    return auroc

def load_weka_filtered_tpairs(fname, min_score=None):
    tested_pairs = [('','',r[0],true_to_1(r[1])) for r in
        ut.load_tab_file(fname)]
    tested_pairs.sort(key=lambda x:x[2], reverse=True)
    if min_score is not None:
        tested_pairs = [t for t in tested_pairs if float(t[2])>=min_score]
    return tested_pairs
    
def preccheck(tested, pchecks=[0.9], total_trues=None):
    recall,precisions = pr(tested)
    recalled = [calc_recall(precisions, p, total_trues) for p in pchecks]
    return recalled

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 2:
        sys.exit("usage: python cv.py filename [total_trues]")
    fname = sys.argv[1]
    pchecks = [0.99,0.90,0.70,0.50]
    recalled = preccheck_wekafiltered(fname, pchecks, None)
    results = '  '.join([(str(r)+'_@_'+str(p)) for (p,r) in
        zip(pchecks,recalled)])
    if nargs == 3:
        total_trues = sys.argv[2]
        recalled = preccheck_wekafiltered(fname, pchecks, total_trues)
        results += '  '.join([(str(r)+'_@_'+str(p)) for (p,r) in
            zip(pchecks,recalled)])
    print fname, results
    
    
    
