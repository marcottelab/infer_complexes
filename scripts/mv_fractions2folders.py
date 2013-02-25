import os
from os.path import abspath
import sys
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut

def maybe_move(fpath, file2folder, remove_final_underscore):
    if not os.path.exists(fpath):
        print "File not found:", fpath
        return
    basename = ut.shortname(fpath)
    if remove_final_underscore:
        basename = ('_'.join(basename.split('_')[:3]) 
                if len(basename.split('_'))>2 else basename)
    if not basename in file2folder:
        print "No mapping for file:", fpath, basename
        return 
    folder = file2folder[basename]
    if not os.path.exists(folder):
        print "Creating directory", folder
        os.mkdir(folder)
    newpath = os.path.join(folder, os.path.split(fpath)[1])
    if os.path.exists(newpath):
        print "File exists:", newpath
    else:
        print "Moving to", newpath
        os.rename(fpath, newpath)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python blah.py files2folders.txt remove_final_underscore{0,1} filename(s)") 
    fname_map = sys.argv[1]
    remove_final_underscore = int(sys.argv[2])
    print "Remove final underscore", "yes" if remove_final_underscore else "no"
    filenames = sys.argv[3:]
    files2folders = dict(ut.load_tab_file(fname_map))
    for f in filenames:
        maybe_move(f, files2folders, remove_final_underscore)
