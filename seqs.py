import operator
import utils as ut
import difflib

def ensembl_prots_to_genes(fname, bar_split=None, second_split=False, 
        only_geneid_on_line=False, pid_replace=False):
    """
    Take a protein sequence file and keep only the longest sequence for each
    gene.  Designed for ensembl fasta sequence downloads.  Purpose is to run
    inparanoid for orthology only on the longest sequence per gene so as to
    have gene-based orthology, which is cleaner to understand.
    bar_split: use 1 for Dd, Sp, leave out for standard ensembl
    for Sp, use second_split=True
    pid_replace only works the first time--don't try it again after replacement
    """
    genes_dict = _longest_seqs_dep(fname, bar_split, 
            only_geneid_on_line=only_geneid_on_line)
    if pid_replace:
        for geneid, lines in genes_dict.items():
            items = lines[0].split(' ')
            protid = items[0].strip().strip('>')
            items[0] = '>' + geneid
            items.append('protein:' + protid)
            lines[0] = ' '.join(items)
    genes_list = reduce(operator.add,[lines for g,lines in genes_dict.items()])
    ut.write_tab_file(genes_list, fname+'_longest')

def prots2genes(fname):
    """
    If there's only one item in the first line, just return a dummy dict
    mapping each id to itself.
    Otherwise, assume the line begins with >GENEID and ends with
    protein:PROTEINID.
    """
    lines = [l for l in ut.load_list(fname) if l[0]=='>']
    if len(lines[0].split())==1:
        return dict([(g,g) for g in [l.strip('>') for l in lines]])
    elif len(lines[0].split(':'))==1:
        # Xl
        return dict([(g,g) for g in [l.split()[0].strip('>') for l in lines]])
    else:
        return dict([(p.split()[-1].split(':')[1], p.split()[0].strip('>'))
                    for p in lines])


def _longest_seqs_dep(fname, bar_split, second_split=False,
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

def load_seqs_from_fasta(fname):
    lol = _load_prots_to_lol(fname)
    return dict([(item[0].split()[0].strip('>'),''.join(item[1:])) 
        for item in lol])

def load_prots_from_fasta(fname):
    """
    Files are in data/sequences/canon.  
    Returns a set since usually I'm searching against it.
    """
    protlines = [l for l in ut.load_list(fname) if l[0]=='>']
    genes = set([l.split(' ')[0].strip('>') for l in protlines])
    return genes

def load_prots_from_fasta_dep(fname):
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

def all_p2g(fs):
    return reduce(ut.dict_quick_merge, [prots2genes(f) for f in fs])

def orth_pid2geneid(fname, p2g):
    lines = ut.load_tab_file(fname)
    def process(lines):
        def replistp2g(pclist):
            return ' '.join([el if i%2 else p2g[el] 
                            for i,el in enumerate(pclist)])
        for n,items in enumerate(lines):
            if n==1:
                yield items
            else:
                newitems = list(items[:2])
                for i in 2,3:
                    newitems.append(replistp2g(items[i].split()))
                yield newitems
    ut.write_tab_file(process(lines), fname+'_fix')

def elut_p2g(fname, p2g, suffix='_fix'):
    lines = ut.load_tab_file(fname)
    def process(lines):
        for items in lines:
            if items[0][0] != '#':
                yield [p2g[items[0]]] + list(items[1:])
            else:
                yield items
    ut.write_tab_file(process(lines), fname+suffix)


class GTrans(object):

    def __init__(self, sp='Hs'):
        lines = ut.load_list_of_lists(ut.proj_path('gene_desc_'+sp))[1:]
        processed = [(l[0],l[1].lower(),l[2] if len(l)>2 else '') for l in lines]
        self.gnames = [(l[1], l[2]) for l in processed]
        self.name2id = dict([(l[1],l[0]) for l in processed])
        self.id2name = dict([(l[0], l[1]) for l in processed])
        self.name2desc = dict(self.gnames)
        self.id2desc = dict([(l[0],l[2]) for l in processed])


    def gfind(self, name):
        name = name.lower()
        return [g for g in self.gnames if g[0].find(name)>-1]

    def find(self, key):
        key=key.lower()
        return [g for g in self.gnames if g[0].find(key)>-1 or
                g[1].find(key)>-1]


def ids2names(ids, gt=None, **kwargs):
    gt = gt or GTrans(**kwargs)
    return [gt.id2name[i] for i in ids]

def names2ids(names, **kwargs):
    gt = GTrans(**kwargs)
    return [gt.name2id[n] for n in names]
