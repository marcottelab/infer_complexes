from __future__ import division
from Struct import Struct

def load_go_ont(fname, filter_namespace='cellular_component'):
    dterms = {}
    term = None
    for r in file(fname,'r'):
        if r.strip()=='[Term]':
            if hasattr(term,'set_is_a'):
                if not filter_namespace or filter_namespace==term.namespace:
                    dterms[term.acc] = (term.name, term.set_is_a)
            term = Struct()
        elif r.strip() and term:
            att = r.split(': ')[0]
            val = ': '.join(r.split(': ')[1:]).strip()
            if att=='id':
                term.acc = val
            elif att=='name':
                term.name = val
            elif att=='namespace':
                term.namespace = val
            elif att=='is_a':
                val = val.split(' ! ')[0]
                if hasattr(term,'set_is_a'):
                    term.set_is_a.add(val)
                else:
                    term.set_is_a = set([val])
    return dterms
