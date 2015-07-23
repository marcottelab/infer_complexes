from __future__ import division
import argparse
import analyze as ana
import ppi_utils as pu
import utils as ut

def holdout_poisson_enrichment(sps_files_dict, cxs_or_ppis, data_out_fname):
    """
    libra: spd = dict([(s,sorted(glob.glob('/home/blakeb/elution/hidden_v35/%s*peps2_FDR001' %s)+glob.glob('/home/blakeb/elution/hidden_v35/%s*peps2_FDR0010' %s))) for s in sps])
    """
    arr = ana.enrichment_array_combined('Hs',sps_files_dict, cxs_or_ppis)
    np.save(data_out_fname, arr)
