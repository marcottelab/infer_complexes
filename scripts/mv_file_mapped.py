import os
from os.path import abspath
import sys
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut

def move(fname, fmap):
    """
    For renaming a file based on a mapping old_fname to new_fname.
    NOT for moving a file to mapped folder.  That's the other script.
    """
    basename = ut.shortname(fname)
    fext = os.path.splitext(fname)[1]
    fdir = os.path.split(fname)[0]
    if basename in fmap:
        newname = os.path.join(fdir,fmap[basename] + fext)
        print "moving", fname, newname
        os.rename(fname, newname)
    else:
        print "not found", fname

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python blah.py mapfile.txt filename(s)") 
    fname_map = sys.argv[1]
    filenames = sys.argv[2:]
    fmap = dict(ut.load_tab_file(fname_map))
    for f in filenames:
        move(f, fmap)
