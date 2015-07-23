from __future__ import division
import argparse
from Bio import SeqIO
from collections import defaultdict
import itertools as it
import numpy as np
import os
import random
import decorators
import orth
import ppi_utils as pu
import seqs
import utils as ut

def load_odict(fname):
    def ogroups_to_odict_list(ogroups):
        orthdict = dict([(p1,[p2 for p2 in og[1]]) for og in ogroups for
                    p1 in og[0]])
        return orthdict
    odict = ogroups_to_odict_list(orth.load_ogroups('','',fname=fname))
    return odict

def load_seqs(fasta_fname):
    records = [x for x in SeqIO.parse(fasta_fname, "fasta")]
    medlen = np.median([len(r.seq) for r in records])
    print "%s: %s sequences, median length %s" % (ut.shortname(fasta_fname),
            len(records), medlen)
    return records

def load_seq_dict(fasta_fname):
    records = load_seqs(fasta_fname)
    sdict = dict([(r.id, r.seq) for r in records])
    return sdict

def all_identities(source_ps, odict_fname, source_fasta, target_fasta,
        target_id_dict_fname=None):
    odict = load_odict(odict_fname)
    if target_id_dict_fname is not None:
        tid_dict = dict([(x[0],x[2]) for x in
            ut.load_lot(target_id_dict_fname)])
        odict = dict([(k,[tid_dict[v] for v in vs]) for k,vs in odict.items()])
    dsource = load_seq_dict(source_fasta)
    dtarget = load_seq_dict(target_fasta)
    pairs = [(s, odict[s][0]) for s in source_ps if s in odict]
    print "%s of %s with orthologs--getting identities" % (len(pairs),
            len(source_ps))
    idents = []
    for s,t in pairs:
        try:
            ident = seqs.percent_identity(dsource[s], dtarget[t])
            idents.append(ident)
        except decorators.TimeoutError, ex:
            print "timeout for %s %s" %(s,t), ex.args
        except Exception, ex:
            print "unknown error for %s %s" % (s,t), ex.args
    print "%s of %s computed successfully" % (len(idents), len(pairs))
    return idents

def ensg_to_ensp_and_park(ppips):
    dhpg = seqs.prots2genes('/Users/blakeweb/Dropbox/complex/data/sequences/canon/Hs.fasta')
    dhgp = ut.dict_inverse(dhpg)
    parkids = ut.load_lol('./orth_similarities/table.Hsapiens/Hsapiens_id.txt')
    ppips_ensp = [dhgp[g] for g in ppips]
    dg2park = dict([(x[2],x[0]) for x in parkids])
    dp2park = dict([(x[1],x[0]) for x in parkids])
    park_ppips_most = [dp2park[p] for p in ppips_ensp if p in dp2park]
    ppips_ensp_rest = [p for p in ppips_ensp if p not in dp2park]
    ppips_ensg_rest = [dhpg[p] for p in ppips_ensp_rest]
    park_ppips_rest = [dg2park[p] for p in ppips_ensg_rest if p in dg2park]
    park_ppips = park_ppips_most + park_ppips_rest
    return park_ppips

def multi_identities(input_fname, out_dir):
    input_list = ut.load_lol(input_fname)
    for desc, prots_fname, source_fasta, odict, target in input_list:
        print "%s, proteins: %s\n source: %s\n odict: %s\ntarget: %s" % (desc,
                prots_fname, source_fasta, odict, target)
        prots = ut.load_list(prots_fname)
        sims = all_identities(prots, odict, source_fasta, target)
        out_fname = os.path.join(out_dir,
                ut.shortname(target).split('.')[0] + "_" + desc + ".txt")
        ut.write_tab_file(sims, out_fname, islist=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fname")
    parser.add_argument("out_dir")
    a = parser.parse_args()
    multi_identities(a.input_fname, a.out_dir)
