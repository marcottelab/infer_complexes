from collections import Counter
import os
from os.path import abspath
import sys
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut

def countfs(fmap, filenames):
    fshorts = [f.split('.')[0] for f in filenames]
    counts = [count_dict_values(fmap, x) for x in fmap.keys(), fshorts]
    print "folder original_counts current_counts"
    output = [(folder, counts[0][folder], counts[1][folder]) for folder in
            sorted(set(fmap.values()))]
    for x in output: print x[0], x[1], x[2]
    print "\n\nfinished:"
    for x in output: 
        if x[1]==x[2]: print x[0]

def count_dict_values(d, keys):
    c = Counter()
    for k in keys:
        if k in d: # some files are renamed, eg 5_8cyto4 to 5_8cyto04. fix.
            c[d[k]] += 1
    return c

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python blah.py mapfile.txt filename(s)") 
    fname_map = sys.argv[1]
    filenames = sys.argv[2:]
    fmap = dict(ut.load_tab_file(fname_map))
    countfs(fmap, filenames)
