import os
from os.path import abspath
import sys
import subprocess
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut

exec_unf = os.path.expanduser(
"~/Dropbox/complex/tools/unfinnigan/perl/Finnigan/bin/uf-mzml -c"
) 
#exec_unf = os.path.expanduser(
#"~/Dropbox/complex/tools/unfinnigan/perl/Finnigan/bin/uf-mzml"
#) # if the first doesn't work, try this (-c is for centroids--depends on how
## the mass spec machine was run
exec_msconvert = os.path.expanduser(
"~/Dropbox/complex/tools/pwiz/msconvert --32 --mzXML"
)

def process(fpath, destdir, do_copysource):
    fname = os.path.split(fpath)[1]
    new_fname = os.path.splitext(fname)[0] + '.mzXML.gz'
    pathout = os.path.join(destdir, new_fname)
    if not os.path.exists(pathout):
        tempdir = pathout + ".tmp"
        if ut.temp_placeholder(tempdir):
            if do_copysource:
                fpath = copy_source(fpath, destdir)
            f_mzml = unf(fpath, destdir)
            if do_copysource:
                os.remove(fpath)
            f_mzxml = msconvert(f_mzml)
            os.remove(f_mzml)
            compress(f_mzxml)
            os.rmdir(tempdir)
        else:
            print "Placeholder exists:", tempdir
    else:
        print "Output exists:", pathout

def copy_source(fpath, destdir):
    ut.run_command('cp %s %s' % (fpath, destdir))
    pathnew = os.path.join(destdir, os.path.split(fpath)[1])
    return pathnew

def unf(pathin, destdir):
    fname = os.path.split(pathin)[1]
    outname = os.path.splitext(fname)[0] + '.mzML'
    outpath = os.path.join(destdir, outname)
    cmd = "%s %s > %s" % (exec_unf, pathin, outpath)
    ut.run_command(cmd)
    return outpath

def msconvert(fin):
    fout = fin.replace('.mzML','.mzXML')
    out_newdir = os.path.splitext(fout)[0]
    cmd = "%s %s -o %s" % (exec_msconvert, fin, out_newdir)
    ut.run_command(cmd)
    match_output = os.path.join(out_newdir, '*.mzXML')
    ut.run_command('mv %s %s' % (match_output, fout))
    os.rmdir(out_newdir)
    return fout

def compress(fname):
    cmd = "gzip %s" % fname
    ut.run_command(cmd)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("usage: python blah.py filename destdir copy_source{0,1}") 
    filename = sys.argv[1]
    destdir = sys.argv[2]
    do_copysource = int(sys.argv[3])
    process(filename, destdir, do_copysource)
