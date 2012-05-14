import random
import sys
import os
import glob
import utils as ut
from Struct import Struct
import complex as co
import elution as el
import fnet
import conf
id_convert_dir = '~/Dropbox/complex/data/sequences/convert'

def shuffled_base(positives, negatives):
    """
    positives, negatives: like [('id1','id2','true'/'false'), ...]
    """
    examples = positives + negatives
    random.shuffle(examples)
    exstruct = Struct(examples=examples, names=['id1','id2','hit'])
    return exstruct

def exstruct_merge_noshuf(exs1, exs2):
    assert exs1.names == exs2.names
    exstruct = Struct(examples=exs1.examples+exs2.examples, names=exs1.names)
    return exstruct

def exstruct_split(exs, nsplit):
    exs1 = Struct(examples=exs.examples[:nsplit], names=exs.names)
    exs2 = Struct(examples=exs.examples[nsplit:], names=exs.names)
    return exs1, exs2

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

def examples_combine_scores(ex_struct, index_start, index_end, reduce_func,
        retain_scores=False, default=0, unknown='?'):
    """
    Combines scores at the specified index columns using reduce_func.
    Put the new score as a new column on the end or replaces with new.
    """
    exs = ex_struct.examples
    def replace_unknowns(lst):
        newlist = []
        for l in lst:
            newval = l if l!=unknown else default
            newlist.append(newval)
        return newlist
    newexs = [e if retain_scores else e[:index_start] + [reduce(reduce_func,
        replace_unknowns(e[index_start:index_end]))] for e in exs]
    names = (ex_struct.names if retain_scores else ex_struct.names[:index_start]) + \
            [ex_struct.names[index_start] + reduce_func.__name__]
    ex_struct.examples = newexs
    ex_struct.names = names

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

def base_examples(key):
    """
    npos=nnegs=None means to use all pos and matching length negs.
    """
    pn_files = ppi_files(key) # order: trainpos, trainneg, testp, testn
    train,test = [shuffled_base(co.pairs(pn_files[p]),co.pairs(pn_files[n]))
                    for (p,n) in [(0,1),(2,3)]]
    ex_struct = exstruct_merge_noshuf(train, test)
    return ex_struct, len(train.examples)

def full_examples(key, elut_fs, scores, species, fnet_gene_dict, suffix='',
                  out_base='weka/', elut_score_cutoff=0.5):
    """
    Key like 'Ce_ensp', 'Hs_uni'. species like 'Hs'.
    Use fnet_gene_dict = -1 to skip functional network.  None means no dict is
        needed. Can supply the dict itself or a string--like 'cep2ceg' or
        'paper_uni2ensg'
    npos=nnegs=None means to use all pos and matching length negs.
    For the set of train_frac, load equal pos and neg.  For remainder (test)
        load nnegs negs.
    """
    # Train and test are merged then split to speed this up 2x
    ex_struct, ntrain = base_examples(key)
    el.score_multi_elfs(ex_struct, elut_fs, scores)
    # Filter out train AND test examples without a score exceeding cutoff
    if elut_fs and elut_score_cutoff is not None:
        ex_struct, ntrain = split_filt_merge(ex_struct, range(3,
                  len(ex_struct.names)), elut_score_cutoff, ntrain)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_struct, species, genedict=fnet_gene_dict)
    out_fname = os.path.join(out_base, species+'_'+key+'_'+suffix+'.arff')
    out_fname = dont_overwrite(out_fname)
    exs_train, exs_test = exstruct_split(ex_struct, ntrain)
    weka_export(exs_train, ut.pre_ext(out_fname,'_train'))
    weka_export(exs_test, ut.pre_ext(out_fname,'_test'))
    print 'ready:', out_fname
    return out_fname

def dont_overwrite(fname):
    if os.path.exists(fname):
        return ut.pre_ext(fname,str(random.randint(0,100)))
    else:
        return fname
    

def predict_all(elut_fs, scores, species, fnet_gene_dict, suffix='',
                  out_base='weka/', elut_score_cutoff=0.5):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    pairs = el.all_filtered_pairs(elut_fs, scores, elut_score_cutoff)
    # examples like [['id1', 'id2', 'true/false'], ...]
    exs = [[p1, p2, '?'] for p1,p2 in pairs]
    ex_struct = Struct(examples=exs,names=['id1','id2','hit'])
    el.score_multi_elfs(ex_struct, elut_fs, scores)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_struct, species, genedict=fnet_gene_dict)
    out_fname = dont_overwrite(os.path.join(out_base,
                                species+'_'+suffix+'.arff'))
    weka_export(ex_struct, out_fname)
    print 'ready:', out_fname
    return out_fname

def split_filt_merge(ex_struct, columns, cutoff, n):
    etrain, etest = exstruct_split(ex_struct, n)
    [filter_scores(exs, columns, cutoff) for exs in [etrain, etest]]
    return exstruct_merge_noshuf(etrain, etest), len(etrain.examples)

def filter_scores(ex_struct, columns, cutoff, missing='?'):
    new_exlist = [e for e in ex_struct.examples if default_max([e[i] for i in columns
        if e[i]!=missing], 0) > cutoff]
    ex_struct.examples = new_exlist

def default_max(numlist, default):
    if numlist==[]: return default
    else: return max(numlist)

def ppi_files(key):
    return [ut.projpath('corum_pairs', key+'_pairs_ppi_'+extra+'.tab') for
        extra in ['train','train_negs','test','test_negs']]

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 5:
        sys.exit("usage: python ml.py train_key elut_file_pattern data_scores species \
            fnet_gene_dict suffix ")
    print sys.argv
    train_key = sys.argv[1]
    elut_files=glob.glob(os.path.expanduser(sys.argv[2]))
    scores = sys.argv[3].split(',')
    species = sys.argv[4]
    fnet_dict = -1 if sys.argv[5]=='-1' else sys.argv[5]
    suffix = sys.argv[6]
    out_path=''
    if train_key == 'full':
        out_fname = predict_all(elut_files, scores, species, fnet_dict, suffix,
            out_path)
    else:
        out_fname = full_examples(train_key, elut_files, scores, species,
            fnet_dict, suffix, out_path)
