import os
import sys

def transpose(d, fin, fout):
    sys.path.append(d+'/..')
    import utils as ut
    #ignore comments, such as last line in spcount output
    endindex = -1 if fin.find('spcount')>-1 else None
    if endindex: print "skipping last line"
    lines = [l for l in ut.load_tab_file(fin)][:endindex]
    cols = ut.zip_exact(*lines) #zip messes up if these files aren't neat
    col2title = cols[1][0].lower()
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
 
