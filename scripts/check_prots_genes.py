from __future__ import division
import os
from os.path import abspath
import sys
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import elution as el
import utils as ut
import seqs

def check(fasta, protq, do_convert):
    p2g = seqs.prots2genes(fasta)
    g2p = ut.dict_inverse(p2g)
    fprots = el.load_elution(protq).prots
    print "checking", ut.shortname(protq)
    print "proteins: %s of %s" % (len([p for p in fprots if p in p2g]),
            len(fprots))
    ngenesfound = len([p for p in fprots if p in g2p])
    print "genes: %s of %s" % (ngenesfound,
            len(fprots))
    if do_convert and ngenesfound < len(fprots):
        print "converting prots to genes:",  protq
        seqs.elut_p2g(protq, p2g)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("usage: python check.py sp.fasta convert{0/1} elut_output ") 
    fasta = sys.argv[1]
    do_convert = int(sys.argv[2])
    elutfs = sys.argv[3:]
    for fname in elutfs:
        check(fasta, fname, do_convert)
