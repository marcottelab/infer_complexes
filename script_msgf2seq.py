from __future__ import division
import sys
import itertools as it
import numpy as np
import utils as ut

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

def msgfbest2sequest_line(r):
    fname_ind_scan_charge = r[0]
    charge = float(r[1])
    msgf_mass = float(r[2])
    msgf_diff = float(r[3])
    peptide = r[4]
    protein = r[5].replace('rv_','rm_')
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
            protein,
            'Z.%s.Z' % peptide
            ]])
    return line_out

#def msgf2sequest_line(r):
    #filename = r[0].replace('.mzXML','')
    #specind, scan = r[1:3]
    #msgf_mass = float(r[4])
    #charge = int(r[6])
    #peptide = r[7].split('.')[1]
    #protein = r[8].replace('rv_','rm_')
    #seq_mass = msgf2seq_mass(msgf_mass, charge)
    #seq_diff = calculate_mass(peptide) - seq_mass
    #line_out = ' '.join([str(i) for i in 
            #[
            #1, # placeholder
            #'.'.join([filename, specind, scan, str(charge)]),
            #seq_mass,
            #'({0:+.05f})'.format(seq_diff),
            #5,5,5,1,"12/345",1234.5,0, # placeholders
            #protein,
            #'Z.%s.Z' % peptide
            #]])
    #return line_out

def msgf2seq_mass(m, z):
    # msgfdb uses observed m/z, sequest uses an estimate of the molecular mass
    # of the observed peptide
    return m * z - z * PROTON

def calculate_mass(peptide):
    return sum([AAS[aa] for aa in peptide]) + WATER

    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: python script_msgf2seq.py filename") 
    filename = sys.argv[1]
    fout = next(it.islice(ut.load_tab_file(filename),2,3))[0].split('.')[0]
    fout += '.mzXML_dta.txt'
    in_gen = ut.load_tab_file(filename)
    in_gen.next(); in_gen.next() # skip 2 lines
    ut.write_tab_file((msgfbest2sequest_line(r) for r in in_gen), fout)
