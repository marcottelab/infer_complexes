def rank_pairs(scores, sample_frac=.1):
    import random
    ranked_scores = []
    for i in range(scores.shape[0]):
        for j in range(scores.shape[1]):
            ranked_scores.append((i,j,scores[i,j]))
    ranked_scores = random.sample(ranked_scores, int(scores.size * sample_frac))
    ranked_scores.sort(key = lambda x: x[2], reverse = True)
    return ranked_scores

def test_pairs(scores, true_pairs, sample_frac=.1):
    # scores: input 2d array of scores for each index-index pair
    # true_pairs: dict: {gene1: set(gene2,gene3,gene4), gene2:
    #   set(gene1,gene5)}
    ranked = rank_scores(scores, sample_frac)
    tested = []
    for i,j,score in ranked:
        hit = int(j in true_pairs.get(i,[]))
        tested.append( (i, j, score, hit) )
    return tested
