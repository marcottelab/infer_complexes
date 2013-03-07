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
    print "proteins: %s of %s" % (len([p for p in fprots if p in p2g]),
            len(fprots))
    ngenesfound = len([p for p in fprots if p in g2p])
    print "genes: %s of %s" % (ngenesfound,
            len(fprots))
    if do_convert and ngenesfound < len(fprots):
        print "converting prots to genes:" protq
        seqs.elut_p2g(protq, p2g)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("usage: python check.py sp.fasta protein_quant_output convert{0/1}") 
    fasta = sys.argv[1]
    protq = sys.argv[2]
    do_convert = int(sys.argv[3])
    check(fasta, protq, do_convert)
