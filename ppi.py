import examples as ex
import complex as co
import seqs
import orth
import utils as ut
import elution as el
from Struct import Struct
import fnet

def learning_examples(species, seqdb, elut_fs, scores=['poisson','wcc','apex'],
                      fnet_file=None, splits=[0,.33,.66,1], neg_ratios=[10,40],
                      ind_cycle=[0,-1], score_cutoff=0.25, base_exstructs=None,
                      pos_splits=None):
    """
    Species: 'Hs', 'Ce', ...
    seqdb: 'uni', 'ensp'
    Use fnet_file = -1 to skip functional network. Use default for it to be
      figured out from species and seqdb as in ut.config.
    Provide base_exstructs as [extr,exte] to use saved examples.
    """
    elut_prots = load_prot_set(elut_fs) # for test set negatives
    if base_exstructs is None:
        ppi_cxs,clean_cxs = load_training_complexes(species, seqdb)
        exstructs = ex.base_examples(ppi_cxs, clean_cxs, elut_prots, splits,
            nratio_train=neg_ratios[0], nratio_test=neg_ratios[1],
            pos_splits=pos_splits, ind_cycle=ind_cycle)
    else:
        exstructs = base_exstructs
    el.score_multi_exs(exstructs, elut_fs, scores, score_cutoff)
    ntest_pos = len([e for e in exstructs[1].examples if e[2]=='true'])
    print exstats(exstructs)
    for exs in exstructs:
        if score_cutoff != -1:
            filter_scores(exs, range(3, len(exs.names)), score_cutoff)
        if fnet_file!=-1:
            fnet.score_examples(exs, species, seqdb, fnet_file)
    print exstats(exstructs)
    return exstructs, ntest_pos #[exstruct_train, exstruct_test]

def load_prot_set(elut_fs):
    return reduce(set.union, (set(el.load_elution(f).prots) for f in elut_fs))

def ex_struct(examples,names):
    return Struct(examples=examples, names=names)

def ex_struct_copy(exs):
    return ex_struct(list(exs.examples),list(exs.names))

def exstats(extr_exte):
    stats = [len([e for e in ex.examples if e[2]==tf]) for ex in extr_exte for tf in ['true','false']]
    return 'train %sP/%sN; test %sP/%sN' % tuple(stats)
    
def load_training_complexes(species, seqdb):
    ppi = co.load_complexes_singleline(ut.proj_path('ppi_cxs'))
    clean = co.load_complexes_multiline(ut.proj_path('clean_cxs'))
    ppi,clean = [convert_complexes(cxs, species, seqdb) for cxs in [ppi,clean]]
    return ppi,clean

def convert_complexes(cxs, species, seqdb):
    #assert species=='Hs' or seqdb=='ensp', "translation not supported"
    if seqdb=='ensp':
        dconv = ut.load_dict_sets(ut.proj_path('convert','Hs_uni2Hs_ensp.tab'))
        ps = seqs.load_prots_from_fasta(ut.proj_path('ensp_fasta','Hs_longest'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    if species!='Hs':
        from_sp = 'Hs_'+seqdb
        to_sp = species+'_'+seqdb
        dconv = orth.odict(from_sp, to_sp)
        ps = seqs.load_prots_from_fasta(ut.proj_path('ensp_fasta',
            species+'_longest'))
        cxs = co.convert_complexes(cxs, dconv, ps)
    return cxs

def filter_scores(ex_struct, columns, cutoff, missing='?'):
    def default_max(numlist, default):
        if numlist==[]: return default
        else: return max(numlist)
    new_exlist = [e for e in ex_struct.examples if default_max([e[i] for i
        in columns if e[i]!=missing], 0) > cutoff]
    ex_struct.examples = new_exlist
    
def predict_all(elut_fs, scores, species, fnet_gene_dict,
                elut_score_cutoff=0.25):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    pairs = el.all_filtered_pairs(elut_fs, scores, elut_score_cutoff)
    print len(pairs), 'total interactions passing cutoff'
    exs = [[p1, p2, '?'] for p1,p2 in pairs]
    ex_struct = Struct(examples=exs,names=['id1','id2','hit'])
    el.score_multi_exs([ex_struct], elut_fs, scores, elut_score_cutoff)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_struct, species, genedict=fnet_gene_dict)
    return ex_struct

def preds_thresh(tested, thresh):
    limit = None
    for i,t in enumerate(tested):
        if t[2]<thresh:
            print 'limit:', i, t
            limit = i
            break
    return limit
    
def cyto_export(tested, fname, negmult=100):
    ut.write_tab_file([(t[0],t[1],myml.rescale(t[2],negmult)) for t in
        tested], fname)
