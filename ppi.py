import examples as ex
import complex as co
import seqs
import orth
import utils as ut
import elution as el
from Struct import Struct
import fnet

def learning_examples(species, seqdb, elut_fs, scores, fnet_gene_dict,
                      splits=[0,.33,.66,1], npos=None, neg_ratios=[10,40],
                      elut_score_cutoff=0.5, exsntrain=None, confirm_redo=True):
    """
    Species: 'Hs', 'Ce', ...
    seqdb: 'uni', 'ensp'
    Use fnet_gene_dict = -1 to skip functional network.  None means no dict is
        needed. Can supply the dict itself or a string--like 'cep2ceg.tab' or
        'paper_uni2ensg.tab'
    npos=nnegs=None means to use all pos and matching length negs.
    For the set of train_frac, load equal pos and neg.  For remainder (test)
        load nnegs negs.
        Provide exsntrain as (exs, ntrain) to use saved examples.
    """
    # Train and test are merged then split to speed this up 2x
    if exsntrain is None:
        ppi_cxs,clean_cxs = load_training_complexes(species, seqdb)
        ex_struct, ntrain = ex.base_examples(ppi_cxs, clean_cxs, splits,
              neg_ratios=neg_ratios, pos_lengths=npos,
                      confirm_redo=confirm_redo)
    else:
        ex_struct, ntrain = exsntrain
    ex_scored = Struct(examples=[e for e in ex_struct.examples],names=[n for n
                      in ex_struct.names])
    el.score_multi_elfs(ex_scored, elut_fs, scores)
    # Filter out train AND test examples without a score exceeding cutoff
    if elut_fs and elut_score_cutoff is not None:
        ex_scored, ntrain = split_filt_merge(ex_scored, range(3,
                  len(ex_scored.names)), elut_score_cutoff, ntrain)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_scored, species, genedict=fnet_gene_dict)
    exs_train, exs_test = exstruct_split(ex_scored, ntrain)
    return exs_train, exs_test

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

def split_filt_merge(ex_struct, columns, cutoff, n):
    def filter_scores(ex_struct, columns, cutoff, missing='?'):
        def default_max(numlist, default):
            if numlist==[]: return default
            else: return max(numlist)
        new_exlist = [e for e in ex_struct.examples if default_max([e[i] for i in
                        columns if e[i]!=missing], 0) > cutoff]
        ex_struct.examples = new_exlist
    def exstruct_merge_noshuf(exs1, exs2):
        assert exs1.names == exs2.names
        exstruct = Struct(examples=exs1.examples+exs2.examples, names=exs1.names)
        return exstruct
    etrain, etest = exstruct_split(ex_struct, n)
    [filter_scores(exs, columns, cutoff) for exs in [etrain, etest]]
    return exstruct_merge_noshuf(etrain, etest), len(etrain.examples)

def exstruct_split(exs, nsplit):
    exs1 = Struct(examples=exs.examples[:nsplit], names=exs.names)
    exs2 = Struct(examples=exs.examples[nsplit:], names=exs.names)
    return exs1, exs2

def predict_all(elut_fs, scores, species, fnet_gene_dict,
                elut_score_cutoff=0.5):
    """
    Same more or less as full_examples above, but produces all predictions in
                  the elution files.
    """
    pairs = el.all_filtered_pairs(elut_fs, scores, elut_score_cutoff)
    # examples like [['id1', 'id2', 'true/false'], ...]
    exs = [[p1, p2, '?'] for p1,p2 in pairs]
    ex_struct = Struct(examples=exs,names=['id1','id2','hit'])
    el.score_multi_elfs(ex_struct, elut_fs, scores)
    if fnet_gene_dict!=-1:
        fnet.score_examples(ex_struct, species, genedict=fnet_gene_dict)
    return ex_struct
