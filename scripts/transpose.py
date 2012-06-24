import os
import sys


def transpose(d, fin, fout):
    sys.path.append(d+'/..')
    import utils as ut
    lines = [l for l in ut.load_tab_file(fin)]
    if lines[-1][0].startswith('#'):
        #ignore comments, such as last line in spcount output
        lines = lines[:-1]
        print "skipping last line"
    cols = ut.zip_exact(*lines) #zip messes up if these files aren't neat
    # _After_ zipping, get rid of the column 1 header--R doesn't like it.
    col0list = list(cols[0])
    print col0list[0][0] 
    assert (col0list[0][0] == '#' or col0list[0] == 'Locus') # make sure we're removing what we should be
    col0list.remove(col0list[0])
    cols[0] = tuple(col0list)
    col2title = cols[1][0].lower()
    # get rid of the total/descr column
    if col2title.find('total') > -1 or col2title.find('descr') > -1:
        cols.remove(cols[1])
        print "removing second column--extraneous"
    ut.write_tab_file(cols, fout)
    
if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 2:
        sys.exit("usage: python transpose.py infile outfile")
    d = os.path.dirname(sys.argv[0])
    fin = sys.argv[1]
    fout = sys.argv[2]
    transpose(d, fin, fout)
 
