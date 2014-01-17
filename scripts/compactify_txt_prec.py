from __future__ import division
import os
import sys
import numpy as np

def compact(d, scoref, precision=3):
    sys.path.append(d+'/..')
    import utils as ut
    compactf = '%s.p%s.txt.gz' % (scoref, precision)
    print compactf, precision
    ascores = np.loadtxt(scoref)
    formatstring = '%' + '0.%se' % precision
    np.savetxt(compactf, ascores, fmt=formatstring, delimiter='\t')

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 2:
        sys.exit("usage: python compactify.py fname precision")
    d = os.path.dirname(sys.argv[0])
    fname = sys.argv[1]
    precision = sys.argv[2]
    compact(d, fname, precision)

