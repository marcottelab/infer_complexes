from __future__ import division

def cv_pairs(scores, true_pairs, genes, sample_frac=.1):
    # scores: input 2d array of scores for each index-index pair
    # true_pairs: dict: {gene1: set(gene2,gene3,gene4), gene2:
    #   set(gene1,gene5)}
    ranked = rank_scores(scores, sample_frac)
    tested = []
    for i,j,score in ranked:
        if i!=j: # remove same-same
            hit = int(genes[j] in true_pairs.get(genes[i],[]))
            tested.append( (i, j, score, hit) )
    return tested

def rank_scores(scores, sample_frac=.1):
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

def auroc(xs,ys):
    auroc = 0
    xprev = 0
    yprev = 0
    for x,y in zip(xs,ys):
        auroc += yprev*(x-xprev)
        xprev = x
        yprev = y
    auroc = auroc / (x*y)
    return auroc
