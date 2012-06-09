import os
import utils as ut

def score_arr(arr, species, seqdb, fnet_file=None, genedict=None):
    # genedict: { exid: set([ensg1, ensg2,...]), exid2: ...}
    if genedict is None:
        dictf = ut.proj_path('convert_net', '%s_%s2%s_net.tab' % (species,
                                                              seqdb, species))
        if os.path.exists(dictf):
            genedict = ut.load_dict_sets(dictf)
    if fnet_file is None:
        fnet_file = ut.config()['fnet_%s_%s' % (species, seqdb)]
    print 'Functional network:', fnet_file, 'Dict:', os.path.basename(dictf), \
        len(genedict) if genedict else 0, 'keys'
    filename = ut.proj_path('fnet_path',fnet_file)
    net = load_net(filename)
    num_items = len(net.items()[0][1])
    names = fnet_names(fnet_file)
    default = ['?']*num_items
    num_hits = 0
    for i,row in enumerate(arr):
        p1,p2 = row['id1'],row['id2']
        scores = scores_pair(p1, p2, net, genedict, default)
        if scores != default:
            num_hits += 1
            for score,name in zip(scores,names):
                if score != '?':
                    row[name] = score
    print num_hits, 'network scores found for ', len(arr), 'pairs'

def fnet_names(fnet_file):
    filename = ut.proj_path('fnet_path',fnet_file)
    return [l[0].strip() if l[0].find('=')==-1 else l[0].split('=')[0].strip()
            for l in ut.load_tab_file(ut.pre_ext(filename,'_names'))]
    
def load_net(filename):
    """
    Output: dict: { ensg1-ensg2: [score1, score2, ...], ensg1-ensg5: ...}
    """
    net = dict([(_idpair(l[0], l[1]), list(l[2:])) for l in
                ut.load_tab_file(filename)])
    return net

def _idpair(id1, id2):
    return '-'.join([id1,id2])

def scores_pair(id1, id2, net, conv2ensg, default):
    """
    Return the maximum score for each score found in corresponding net
    """
    # Our conversion dict is sets: get a list for each id in the pair
    id1s = conv2ensg.get(id1,[]) if conv2ensg else [id1]
    id2s = conv2ensg.get(id2,[]) if conv2ensg else [id2]
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

def score_examples(ex_struct, species, seqdb, fnet_file=None, genedict=None):
    # exs: struct with exs.names: ['score1', 'score2', ...]
    # and exs.examples: [[id1, id2, 'true/false', score1, score2, ...], ..]
    # genedict: { exid: set([ensg1, ensg2,...]), exid2: ...}
    if genedict is None:
        dictf = ut.proj_path('convert_net', '%s_%s2%s_net.tab' % (species,
                                                              seqdb, species))
        if os.path.exists(dictf):
            genedict = ut.load_dict_sets(dictf)
    if fnet_file is None:
        fnet_file = ut.config()['fnet_%s_%s' % (species, seqdb)]
    print 'Functional network:', fnet_file, 'Dict:', os.path.basename(dictf), \
        len(genedict) if genedict else 0, 'keys'
    filename = ut.proj_path('fnet_path',fnet_file)
    net = load_net(filename)
    num_items = len(net.items()[0][1])
    out_examples = []
    default = ['?']*num_items
    num_hits = 0
    for ex in ex_struct.examples:
        scores = scores_pair(ex[0], ex[1], net, genedict, default)
        out_examples.append(ex + scores)
        if scores != default: num_hits += 1
    print num_hits, 'network scores found for ', len(out_examples), 'pairs'
    ex_struct.examples = out_examples
    ex_struct.names += [l[0] for l in
                        ut.load_list_of_lists(ut.pre_ext(filename,'_names'))]
    
def munge_original(fdata, column_inds, fnames, fout, first_names=1):
    """
    Keep selected columns, replace 'NA' with '?', remove empty rows.
    Do not include 0 or 1 for ids--they are kept automatically.
    For column inds, start with 0 for scores.
    Keep the same columns from the fnames file so I have a record of it.
    """
    out = []
    default = ['?'] * len(column_inds)
    for l in ut.load_tab_file(fdata):
        ids = list(l[:2])
        newdata = [l[i+2] if l[i+2]!='NA' else '?' for i in range(len(l)) if i
            in column_inds]
        if newdata != default:
            out.append(ids + newdata)
    ut.write_tab_file(out, fout)
    names = [l for i,l in enumerate( list( ut.load_tab_file(
        fnames))[first_names:]) if i in column_inds]
    ut.write_tab_file(names, ut.pre_ext(fout, '_names')) 
