import operator
import utils as ut
import difflib

def ensembl_prots_to_genes(fname, bar_split=None, second_split=False, only_geneid_on_line=False):
    """
    Take a protein sequence file and keep only the longest sequence for each
    gene.  Designed for ensembl fasta sequence downloads.  Purpose is to run
    inparanoid for orthology only on the longest sequence per gene so as to
    have gene-based orthology, which is cleaner to understand.
    bar_split: use 1 for Dd, Sp, leave out for standard ensembl
    for Sp, use second_split=True
    """
    genes_dict = _longest_seqs(fname, bar_split, only_geneid_on_line=only_geneid_on_line)
    genes_list = reduce(operator.add,[lines for g,lines in genes_dict.items()])
    ut.write_tab_file(genes_list, fname+'_longest')

def _longest_seqs(fname, bar_split, second_split=False,
                 only_geneid_on_line=False):
    prots = _load_prots_to_lol(fname)
    genes_dict = {}
    for p in prots:
        if only_geneid_on_line:
            geneid = p[0]
        else:
            if bar_split:
                geneid = p[0].split('|')[bar_split].strip()
                if second_split:
                    geneid = geneid.split('.')[0]
            else:
                geneid = p[0].split('gene:')[1].split('transcript:')[0].strip()
        seq_length = sum([len(pi) for pi in p[1:]])
        if seq_length > genes_dict.get(geneid,(0,''))[0]:
            genes_dict[geneid] = (seq_length,p)
    # get rid of the length before returning
    genes_dict = dict([(g,p[1]) for g,p in genes_dict.items()])
    return genes_dict

def _load_prots_to_lol(fname):
    prots = ut.load_list(fname)
    prots_clean = []
    for line in prots:
        if line[0] == '>': 
            prots_clean.append([line])
        else:
            prots_clean[-1].append(line)
    return prots_clean

def load_prots_from_fasta(fname):
    """
    Files are in data/sequences/canon.  All so far can be split by both space
    and |.
    Returns a set since usually I'm searching against it.
    """
    protlines = [l[1:] for l in ut.load_list(fname) if l[0]=='>']
    prots = set([l.split(' ')[0].split('|')[0] for l in protlines])
    return prots

def cuihong_fasta_to_clean(fname, outname):
    """
    Get rid of all the reverse or shuffleds ('rm' instead of 'sp') and anything
    else that doesn't start with '>sp'.  Keep only the uniprot identifier.
    """
    lol = _load_prots_to_lol(fname)
    good_seqs = [['>'+p[0].split('|')[1]]+p[1:] for p in lol 
            if p[0][:3] == '>sp' or p[0][:3]== '>tr']
    ut.write_tab_file([i for l in good_seqs for i in l ], outname) #flatten

def similarities(dseqs):
    """
    Input: dseqs = { geneid1: 'MPOIHUK...', }
    """
    allsims = [];
    for i, (aid,aseq) in enumerate(dseqs.items()):
        sims = []
        for bid,bseq in dseqs.items():
            diff = difflib.SequenceMatcher(a=aseq,b=bseq)
            if aid!=bid and diff.real_quick_ratio() > .90 \
                    and diff.quick_ratio() > .90:
                sims.append((bid,diff.ratio()))
        print i, len(sims), max([s[1] for s in sims]) if sims else 0
        allsims.append((aid,sims))
    return allsims

def ensp_clean_chroms(protlol, set_badchs, start=0, end=0):
    end = len(list(set_badchs)[0]) if start==0 and end==0 else end
    protfix = [plines for plines in protlol
        if plines[0].split(' ')[2].split(':')[2][start:end] not in set_badchs]
    return protfix
