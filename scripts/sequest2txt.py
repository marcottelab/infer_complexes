import os
from os.path import abspath
import subprocess
import random
import sys
from collections import defaultdict
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut

def add_header(seq_header, fname):
    #fout = os.path.join(usedir, fout + '.mzXML_dta.txt')
    fout = fname.split('.')[0] + '.mzXML_dta.txt'
    subprocess.call('cat %s %s > %s' % (abspath(seq_header),
        abspath(fname), abspath(fout)), shell=True)
    os.remove(fname)
    return fout

def keep_unique_lines_by_column(fnames, column=0):
    lines_dict = collect_dict((line for f in fnames 
        for line in ut.load_tab_file(f)))
    return (random.choice(list(lines)) for key,lines in lines_dict.items())

def collect_fnames(filenames):
    return collect_dict(filenames, 
            selector=lambda x:x.split('/')[-1].split('.')[0])

def collect_dict(items, selector=lambda x:x[0]):
    d = defaultdict(set)
    for item in items:
        d[selector(item)].add(item)
    return d

def write_combined(fnames):
    output = keep_unique_lines_by_column(fnames)
    fout = '.'.join(fnames[0].split('.')[:-2]) + '.combined'
    ut.write_tab_file(output,fout)
    return fout

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python script_sequest2txt.py seq_header filename(s)") 
    seq_header = sys.argv[1]
    filenames = sys.argv[2:]
    dict_fnames = collect_fnames(filenames)
    for key,fnames in dict_fnames.items():
        fout = write_combined(list(fnames))
        print "Combined to", fout
        fout = add_header(seq_header, fout)
        print "Added header to", fout
        fout_name = os.path.split(fout)[1]
        newdir = os.path.splitext(fout)[0] 
        assert os.path.exists(newdir), "Dest dir %s does not exist." % newdir
        # move to appropriate folder
        dest = os.path.join(newdir,fout_name)
        #print 'moving:', fout, dest
        os.rename(fout, dest)
