from seqs import GTrans
import utils as ut

def all_ev(gnames, arr_ev, scores):
    """
    Given a list of gene names, return the rows of evidence for any pairs.
    Filter only to those pairs present in scores (weighted predicted
    interaction list--probably the cx_filtered version), if provided.
    """
    gt = GTrans()
    ids = set([gt.name2id[n] for n in gnames])
    if scores:
        pairs = dict([((s[0],s[1]),s[2]) for s in scores 
            if s[0] in ids and s[1] in ids])
        return ([tuple([gt.id2name[i] for i in r[0],r[1]] +
            [pairs[(r[0],r[1])]] + list(r)[3:]) for r in arr_ev if (r[0],r[1])
            in pairs])
    else: 
        return [tuple([gt.id2name[i] for i in r[0],r[1]]+list(r)[3:]) for r in arr_ev if r[0] in ids and r[1] in ids]

def ev_output(gnames, arr_ev, scores=None):
    ev_lot = all_ev(gnames, arr_ev, scores)
    labels = ['gene1','gene2','total_score'] + list(arr_ev.dtype.names[3:])
    return zip(labels, *ev_lot)

