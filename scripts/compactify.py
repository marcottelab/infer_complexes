import os
import sys
import numpy as np

def compact(d, scoref, dtype='f2'):
    sys.path.append(d+'/..')
    import utils as ut
    compactf = '%s.%s.pyd' % (scoref, dtype)
    print compactf, dtype
    ascores = np.loadtxt(scoref, dtype)
    ut.savepy(ascores, compactf)

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 2:
        sys.exit("usage: python compactify.py fname dtype")
    d = os.path.dirname(sys.argv[0])
    fname = sys.argv[1]
    dtype = sys.argv[2]
    compact(d, fname, dtype)

