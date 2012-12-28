from __future__ import division
import itertools as it
import numpy as np
import subprocess
from os.path import abspath
import os
import sys
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

def msgfbest2sequest_line(r, prots2genes):
    fname_ind_scan_charge = r[0]
    charge = float(r[1])
    msgf_mass = float(r[2])
    msgf_diff = float(r[3])
    peptide = r[4]
    protein = r[5]
    protid = protein[3:] if protein.startswith(DECOY) else protein
    geneid = prots2genes[protid]
    seq_protein = SEQ_DECOY + geneid if protein.startswith(DECOY) else geneid
    seq_mass = msgf2seq_mass(msgf_mass, charge)
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

def msgf2seq_mass(m, z):
    # msgfdb uses observed m/z, sequest uses an estimate of the molecular mass
    # of the observed peptide
    return m * z - z * PROTON

def calculate_mass(peptide):
    return sum([AAS[aa] for aa in peptide]) + WATER

    
def msgf2seq_file(filename, fasta_file, seq_header):
    # Get the sample filename from the first item of the third line
    usedir = os.path.split(filename)[0] 
    fout = next(it.islice(ut.load_tab_file(filename),2,3))[0].split('.')[0]
    fout = os.path.join(usedir, fout + '.mzXML_dta.txt')
    in_gen = ut.load_tab_file(filename)
    in_gen.next(); in_gen.next() # skip 2 lines
    p2g = seqs.prots2genes(fasta_file)
    temp_fout = fout + '_tmp'
    ut.write_tab_file((msgfbest2sequest_line(r,p2g) for r in in_gen),
            temp_fout)
    subprocess.call('cat %s %s > %s' % (abspath(seq_header),
        abspath(temp_fout), abspath(fout)), shell=True)
    os.remove(temp_fout)
    return fout

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python script_msgf2seq.py fasta_file seq_header filename(s)") 
    fasta_file = sys.argv[1]
    seq_header = sys.argv[2]
    filenames = sys.argv[3:]
    for f in filenames:
        fout = msgf2seq_file(f, fasta_file, seq_header)
        # move to appropriate folder
        currdir, fout_name = os.path.split(fout)
        new_dirname = fout_name.split('.')[0] + '.mzXML_dta'
        dest = os.path.join(currdir,new_dirname,fout_name)
        #print 'moving:', fout, dest
        os.rename(fout, dest)

