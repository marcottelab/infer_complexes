import os
from os.path import abspath
import sys
import shutil
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import elution as el
import utils as ut

MSB_EXT = '.prot_count_uniqpeps_FDR0010'
PQ_FILE = 'ISO_QUAN/quantification-protein-intensity.txt'
PQ_NEW = '_pqmsb.tab'
PQ_FILT = '_pqmsb_filtmsb.tab'

def process(proj_dir, msb_out_dir, dirnames):
    """
    If a single dirname, just process.
    If multiple, merge then process.
    """
    pq_new_path = os.path.join(proj_dir, proj_name+PQ_NEW)
    if len(dirnames) > 1:
        merge(proj_dir, dirnames, pq_new_path)
    else:
        pq_old_path = os.path.join(proj_dir, PQ_FILE)
        shutil.copy(pq_old_path, pq_new_path)
    pq_filt_path = msb_filter(proj_dir, msb_out_dir, pq_new_path)

def msb_filter(proj_dir, msb_out_dir, pq_new_path):
    """
    Filter the pepquant output by keeping only values with spectral counts in
    the msblender output.
    """
    proj_name = ut.shortname(proj_dir)
    msb_quant_file = os.path.join(msb_out_dir, proj_name+MSB_EXT)
    pq_elut, msb_elut = [el.load_elution(f) for f in pq_new_path,
            msb_quant_file]
    pq_filt = el.filter_matching_elution(pq_elut, msb_elut)
    pq_filt_path = pq_new_path.replace(PQ_NEW, PQ_FILT)
    el.write_elution(pq_filt, pq_filt_path)
    return pq_filt_path
    

def merge(proj_dir, dirnames, pq_new_path):
    """
    Combine pepquant quantitation from project_1 (etc) PQ_FILE into
    project+PQ_NEW.
    """
    if not os.path.exists(proj_dir):
        os.mkdir(proj_dir)
    proj_name = ut.shortname(proj_dir)
    assert not os.path.exists(pq_new_path), "%s exists. Exiting." % pq_new_path
    dirnames = sort_numbered(dirnames)
    print "Sorted dirnames:", dirnames
    pq_files = [os.path.join(d,PQ_FILE) for d in dirnames]
    eluts = (el.load_elution(f) for f in pq_files)
    merged = reduce(el.combine_elutions, eluts)
    el.write_elution(merged, pq_new_path)

def sort_numbered(filenames):
    pairs = [(f, int(f.split('_')[-1])) for f in filenames]
    pairs.sort(key=lambda x: x[1])
    return pairs

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python pepquant_post.py proj_dir msb_out_dir directory(s)") 
    proj_dir = sys.argv[1]
    msb_out_dir = sys.argv[2]
    dirnames = sys.argv[3:]
    print "Directories:", dirnames
    process(proj_dir, msb_out_dir, dirnames)
