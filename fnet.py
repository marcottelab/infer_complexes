import os
import utils as ut
import orth
import pairdict as pd
import itertools as it

def score_arr(arr, species, ext_key):
    ext_file = ut.config()[ext_key]
    print 'External Data:', ext_file
    genedict = convdict_from_fname(species, ext_file)
    filename = ut.proj_path('fnet_path', ext_file)
    net = load_net(filename)
    print 'Ext data size:', len(net)
    nitems = len(net.items()[0][1])
    names = fnet_names(ext_file) if nitems>1 else [ext_key]
    default = ['?']*nitems
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

def convdict_from_fname(species, ext_file):
    # Doesn't yet work for the general case of possibly needing to go two
    # steps--to the new species, then to a new seqdb
    totype = '_'.join(ext_file.split('/')[-1].split('_')[:2]) #Hs_entrez;Dm_fbgn
    # If there's no matching conversion file, assume it's not needed.
    genedict = None
    try:
        genedict = orth.convert_dict(species, totype)
        print 'Conversion file:', species, totype, len(genedict), 'keys'
    except IOError as e:
        print 'No external conversion file:', species, totype, e.strerror
    return genedict

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
    nscores = len(scores)
    if nscores==0:
        return default
    elif nscores==1:
        return scores[0]
    else:
        return max_scores(scores)

def max_scores(scores):
    # TODO: fix this after initial testing
    print "fix max scores--this one was length ", len(scores)
    return scores[0]

def munge_original(fdata, column_inds, fnames, fout, first_names=1):
    """
    Keep selected columns, replace 'NA' with '?', remove empty rows.
    Ids (first 2 columns) are kept automatically.
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

def munge_malov(fdata):
    # load from proper columns
    cxs = {}
    for line in ut.load_tab_file(fdata):
        g,c = line[:2]
        cxs.setdefault(c,set([])).add(g)
    # remove (many) singletons
    for c,gset in cxs.items():
        if len(gset) < 2: 
            del cxs[c]
    ints = pd.PairDict([])
    # interpret "approved"/"provisional"/"temporary"
    def scorec(c):
        if c[0] == 'A':
            return 10
        elif c[0] == 'P':
            return 3
        elif c[0] == 'T':
            return 1
        else:
            print c[0]
            return 1
    for c,gset in cxs.items():
        score = scorec(c)
        for pair in it.combinations(gset,2):
            assert not ints.contains(pair), "ints contains %s" % pair[0]+pair[1]
            ints.append(pair, score)
    return ints
