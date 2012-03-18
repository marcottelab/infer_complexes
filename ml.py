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
    while len(negatives) < K:
        a = random.choice(items) 
        b = random.choice(items)
        newpair = (a,b)
        if newpair not in truepairs and (b,a) not in truepairs:
            negatives.append(newpair)
    return negatives
