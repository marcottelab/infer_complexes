from collections import defaultdict
import utils as ut

def make_prot_counts(f_peplist, only_uniques=False):
    peplist = ut.load_list_of_lists(f_peplist)
    exclude_peps = non_unique_peps(peplist) if only_uniques else set([])
    prot_counts = prot_counts(peplist, exclude_peps)
    samples = sorted(list(set(ut.i1(peplist))))
    suffix = 'prot_count' + ('_uniqpeps' if only_uniques else '')
    f_counts = f_peplist.replace('pep_list',suffix)


def prot_counts(peplist, exclude_peps=set([])):
    """
    Exclude_peps: set of peptides to exclude, probably from non_unique_peps.
    """
    dprots = defaultdict(lambda: defaultdict(ut.constant_factory(0.00)))
    for prot, sample, pep, count in peplist:
        if pep not in exclude_peps:
            dprots[prot][sample] += float(count)
    return dprots

def non_unique_peps(peplist):
    """
    Peplist: list-of-lists read in from an msblender pep_list file
    """
    pep2prots = defaultdict(set)
    for prot,sample,pep,count in peplist:
        pep2prots[pep].add(prot)
    non_unique = set([pep for pep,prots in pep2prots.items() if len(prots)>1])
    return non_unique

