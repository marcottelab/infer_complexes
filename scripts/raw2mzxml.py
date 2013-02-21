import os
from os.path import abspath
import sys
import subprocess
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut

exec_unf = os.path.expanduser(
"~/Dropbox/complex/tools/unfinnigan-14107a88926c/ perl/Finnigan/bin/uf-mzml -c"
)
exec_msconvert = os.path.expanduser(
"""
~/Dropbox/complex/tools/pwiz-bin-linux-x86_64-gcc42-release-3_0_4268/msconvert 
--32 --mzXML
"""
)

def process(fpath, destdir):
    f_mzml = maybe_unf(fpath, destdir)
    maybe_msconvert(f_mzml)
    f_mzxml = os.remove(f_mzml)
    compress(f_mzxml)

def maybe_unf(pathin, destdir):
    fname = os.path.split(pathin)[1]
    outpath = os.path.join(destdir, fname)
    cmd = "%s %s > %s" % (exec_unf, pathin, pathout)
    check_run(cmd, outpath)
    return fout

def check_run(cmd, pathout):
    if not os.path.exists(pathout):
        tempdir = pathout + ".tmp"
        if ut.temp_placeholder(tempdir):
            print cmd
            subprocess.call(cmd)
            os.rmdir(tempdir)

def maybe_msconvert(fin):
    fout = fin.replace('.mzML','.mzXML')
    cmd = "%s %s" % (exec_msconvert, fin)
    check_run(cmd, fout)
    return fout

def compress

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python blah.py filename destdir") 
    filename = sys.argv[1]
    dest = sys.argv[2]
    process(filename, dest)
