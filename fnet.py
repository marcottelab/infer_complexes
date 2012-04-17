import os
import utils as ut

def score_examples(ex_struct, species='Hs', genedict=None):
    # exs: struct with exs.names: ['score1', 'score2', ...]
    # and exs.examples: [[id1, id2, 'true/false', score1, score2, ...], ..]
    # genedict: { exid: set([ensg1, ensg2,...]), exid2: ...}
    net = load_net(species)
    num_items = len(net.items()[0][1])
    out_examples = []
    for ex in ex_struct.examples:
        out_examples.append(ex + scores_pair(ex[0], ex[1], net, genedict,
                                             num_items))
    ex_struct.examples = out_examples
    ex_struct.names += [l for l in
                       ut.load_tab_file(filename(species,which='names'))][0][2:]
    return ex_struct

def load_net(species):
    """
    Output: dict: { ensg1-ensg2: [score1, score2, ...], ensg1-ensg5: ...}
    """
    lines = ut.load_list_of_lists(filename(species))
    net = dict([(_idpair(l[0], l[1]), list(l[2:])) for l in lines])
    return net


def _idpair(id1, id2):
    return '-'.join([id1,id2])

def scores_pair(id1, id2, net, conv2ensg, num_items):
    """
    Return the maximum score for each score found in corresponding net
    """
    default = ['?']*num_items
    # Our conversion dict is sets: get a list for each id in the pair
    id1s = conv2ensg.get(id1,[]) if conv2ensg else []
    id2s = conv2ensg.get(id2,[]) if conv2ensg else []
    pairs = [(a,b) for a in id1s for b in id2s] + [(b,a) for a in id1s for b in id2s]
    if len(pairs)==0:
        return default
    scores = [net[_idpair(p1,p2)] for p1,p2 in pairs if _idpair(p1,p2) in net]
    if len(scores)==0:
        return default
    elif len(scores)==1:
        return scores[0]
    else:
        return max_scores(scores)

def max_scores(scores):
    # TODO: fix this after initial testing
    print "fix max scores--this one was length ", len(scores)
    return scores[0]

#def net2entrez(net_id):
    """
    Functional net ids are just the nonzero end portions of ensemble gene ids.
    Ens gene ids are ENSG and 11 numbers.
    Turns out this was wrong.  They're entrez ids, not ensembl.
    """
#    return 'ENSG'+str(net_id).zfill(11)
        

def filename(species, which='data'):
    base = os.path.expanduser('~/Dropbox/complex/data/functional/')
    if which == 'data':
        return os.path.join(base, species, species+'_filtered.tab')
    elif which == 'names':
        return os.path.join(base, species, 'selected_cols.txt')
        

