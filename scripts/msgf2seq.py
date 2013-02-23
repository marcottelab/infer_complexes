from __future__ import division
import itertools as it
import numpy as np
import subprocess
from os.path import abspath
import os
import sys
sys.path.append(os.path.dirname(abspath(__file__))+'/../')
import utils as ut
import seqs

PROTON = 1.00727646677
WATER = 18.01048 # mass of a water due to no peptide bonds on the end

AAS = dict(
    A = 71.03711,
    R = 156.10111,
    N = 114.04293,
    D = 115.02694,
    C = 160.030654, # iodoacetamide treatment
    #C = 103.00919,
    E = 129.04259,
    Q = 128.05858,
    G = 57.02146,
    H = 137.05891,
    I = 113.08406,
    L = 113.08406,
    K = 128.09496,
    M = 131.04049,
    F = 147.06841,
    P = 97.05276,
    S = 87.03203,
    T = 101.04768,
    W = 186.07931,
    Y = 163.06333,
    V = 99.06841
    )

DECOY = 'rv_'
SEQ_DECOY = 'rm_'

searches = {
        'logSpecProb_hit_list_best': 'msgfdb',
        'MQscore_hit_list_best': 'inspect', 
        'xcorr_hit_list_best': 'tide' 
        }

def msgfbest2sequest_line(r, prots2genes, search):
    fname_ind_scan_charge = r[0]
    charge = float(r[1])
    msgf_mass = float(r[2])
    msgf_diff = float(r[3])
    peptide = r[4]
    protein = r[5]
    protid = protein[3:] if protein.startswith(DECOY) else protein
    geneid = prots2genes[protid] if prots2genes is not None else protid
    seq_protein = SEQ_DECOY + geneid if protein.startswith(DECOY) else geneid
    seq_mass = msgf2seq_mass(msgf_mass, charge, search)
    seq_diff = calculate_mass(peptide) - seq_mass
    seq_diff_temp = seq_diff
    for i in range(5):
        if abs(seq_diff_temp) < .15:
            seq_diff = seq_diff_temp
            break
        else:
            seq_diff_temp += -np.sign(seq_diff_temp)*PROTON
    line_out = ' '.join([str(i) for i in 
            [
            1, # placeholder
            fname_ind_scan_charge,
            seq_mass,
            '({0:+.05f})'.format(seq_diff),
            5,5,5,1,"12/345",1234.5,0, # placeholders
            seq_protein,
            'Z.%s.Z' % peptide
            ]])
    return line_out

def msgf2seq_mass(m, z, search):
    # msgfdb and inspect use observed m/z. 
    # tide uses something weird.
    # sequest uses an estimate of the molecular mass of the observed peptide
    if search in ['msgfdb','inspect']:
        return m * z - z * PROTON
    elif search == 'tide':
        return m * z

def calculate_mass(peptide):
    return sum([AAS[aa] for aa in peptide]) + WATER

    
def msgf2seq_file(filepath, fasta_file, msb_psms):
    """
    msb_psms: set of spectid_peptidesequence
    """
    def parse_spec_pep_row(r):
        # get spec_pep from _best file format
        parsed = '_'.join(r[0].split('.')[:2] + [r[4]])
        #print parsed
        return parsed
    usedir,fin = os.path.split(filepath)
    # Get the sample filename from the first item of the third line
    fout = next(it.islice(ut.load_tab_file(filepath),2,3))[0].split('.')[0]
    in_gen = ut.load_tab_file(filepath)
    in_gen.next(); in_gen.next() # skip 2 lines
    p2g = seqs.prots2genes(fasta_file)
    fout = os.path.join(usedir, '.'.join([fout, fin.split('.')[-1] ,
        'sequestformat']))
    search = searches[filepath.split('.')[-1]]
    print "Converting/filtering; Search:", search
    output = (msgfbest2sequest_line(r,p2g, search) for r in in_gen 
            if parse_spec_pep_row(r) in msb_psms)
    print "Writing", fout
    ut.write_tab_file(output, fout)
    return fout

def parse_msb_psms(fname):
    item1s = (line[0] for line in ut.load_tab_file(fname))
    # ex: WAN110811_HCW_HEK293NE_P1D08.01387.2.SGNLTEDDKHNNAK
    item1s.next() # skip 1 line
    spect_pep = ('_'.join([sample,spect,pep]) for sample,spect,_,pep in 
            (i1.split('.') for i1 in item1s))
    return set(spect_pep)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python script_msgf2seq.py fasta_file msb_psm_file filename(s)") 
    fasta_file = sys.argv[1]
    msb_psm_file = sys.argv[2]
    filenames = sys.argv[3:]
    print "Loading msblender output", msb_psm_file
    msb_psms = parse_msb_psms(msb_psm_file)
    #print "msb psms 0:", list(msb_psms)[0]
    for f in filenames:
        print "Loading search output", f
        fout = msgf2seq_file(f, fasta_file, msb_psms)

