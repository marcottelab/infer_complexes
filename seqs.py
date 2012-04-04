import operator
import utils as ut

    

def ensembl_prots_to_genes(fname, bar_split=None, second_split=False, only_geneid_on_line=False):
    """
    Take a protein sequence file and keep only the longest sequence for each
    gene.  Designed for ensembl fasta sequence downloads.  Purpose is to run
    inparanoid for orthology only on the longest sequence per gene so as to
    have gene-based orthology, which is cleaner to understand.
    bar_split: use 1 for Dd, Sp, leave out for standard ensembl
    for Sp, use second_split=True
    """
    def longest_seqs(fname, bar_split, second_split=False,
                     only_geneid_on_line=False):
        def load_prots(fname):
            prots = ut.load_list(fname)
            prots_clean = []
            for line in prots:
                if line[0] == '>': 
                    prots_clean.append([line])
                else:
                    prots_clean[-1].append(line)
            return prots_clean
        prots = load_prots(fname)
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
        return genes_dict
    genes_dict = longest_seqs(fname, bar_split, only_geneid_on_line=only_geneid_on_line)
    genes_list = reduce(operator.add,[g[1][1] for g in genes_dict.items()])
    ut.write_tab_file(genes_list, fname+'_longest')
    
def load_prots_from_fasta(fname):
    """
    Files are in data/sequences/canon.  All so far can be split by both space
    and |.
    """
    protlines = [l[1:] for l in ut.load_list(fname) if l[0]=='>']
    prots = [l.split(' ')[0].split('|')[0] for l in protlines]
    return prots
    
    
