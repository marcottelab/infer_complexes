import random
import utils as ut

def example_from_interaction_matrix(trueints, mat, labels, default='?'):
    # trueints: interacting partners: [(a,b),(a,d),(b,c),...]
    # mat: N by N matrix of interaction scores
    # labels: length N list of labels for rows/columns (same)
    labeldict = ut.list_inv_to_dict(labels)
    scores = []
    for a,b in trueints:
        if a in labels and b in labels:
            scores.append(mat[labeldict[a],labeldict[b]])
        else:
            scores.append(default)
    return scores

def negative_pairs(items, K, truepairs):
    negatives = []
    dict_pos = full_dict_from_pairs(truepairs)
    while len(negatives) < K:
        a = random.choice(items) 
        b = random.choice(items)
        newpair = (a,b)
        if a not in dict_pos[b] and b not in dict_pos[a]:
            negatives.append(newpair)
    return negatives

def full_dict_from_pairs(pairs):
    d1 = ut.dict_sets_from_tuples(pairs)
    d2 = ut.dict_sets_from_tuples([(p[1],p[0]) for p in pairs])
    return ut.dict_combine_sets(d1,d2)
