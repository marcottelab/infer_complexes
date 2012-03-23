import random
import utils as ut

def interaction_examples(positives, mats_labels, Kneg=None):
    """
    Generate a list of [id1, id2, hit, score1, score2, ...]
    positives: list of true interaction pairs
    mats_labels: list of tuples of [(score_matrix, row_labels), ..]
    Kneg: number of negative examples to generate from items in the pairs
    example:
    ex1 =
    ml.interaction_examples(random.sample(corumpairsen,1000),[(ellong.corr,
    ellong.prots), (elall.corr, elall.prots)], 1000)
    """
    examples = []
    Kneg = Kneg if Kneg is not None else len(positives)
    negatives = negative_pairs(positives, Kneg)
    def unzip(pairlist):
        return [p[0] for p in pairlist], [p[1] for p in pairlist]
    for (pairs,iscomplex) in ((positives, ['true']*len(positives)), (negatives, ['false']*Kneg)):
        scores = [scores_interaction_matrix(pairs, ml[0], ml[1]) for ml in
                  mats_labels]
        p0,p1 = unzip(pairs)
        examples += zip(p0, p1, iscomplex, *scores)
        #examples += zip(iscomplex, *scores)
    random.shuffle(examples)
    return examples

def weka_export(examples, filename, startindex=3):
    # skip the first 2 items specifiying the interaction
    # Note that I"m taking iscomplex(true/false) and putting it at the end
    # per weka convention
    f = open(filename,'w')
    f.write('@RELATION complexes\n\n')
    for i in range(len(examples[0])-3):
        f.write('@ATTRIBUTE score'+str(i)+'\treal\n')
    f.write('@ATTRIBUTE iscomplex\t{true,false}\n')
    f.write('\n@DATA\n\n')
    for ex in examples:
        f.write(', '.join([str(val) for val in ex[3:]]+[ex[2]])+'\n')
    f.write('%\n%\n%\n')
    f.close()

def scores_interaction_matrix(pairs, mat, labels, default='?'):
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

def negative_pairs(truepairs, K):
    # TODO: speed up duplicate checking
    items = list(reduce(set.union,[set(p) for p in truepairs]))
    negatives = []
    def pair_in_list(pair, pairlist):
        return pair in pairlist or (pair[1],pair[0]) in pairlist
    while len(negatives) < K:
        a = random.choice(items) 
        b = random.choice(items)
        newpair = (a,b)
        if a!=b and not pair_in_list(newpair,truepairs) and not \
           pair_in_list(newpair, negatives):
            negatives.append(newpair)
    return negatives
