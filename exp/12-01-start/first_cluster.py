import ml_complex as ml
import plotting as pl
import load_elution as le

# elut_files = [
#     ('Hs_paper','data/12-01-start/elut_hs.tab'),
#     ('Dm_0','data/12-01-start/elut_dm.tab')
#     ('Dm_1','data/12-01-start/elut_dm.tab')
#     ]

elut_files = [
     ('Hs0', '../../data/12-01-start/SupplementaryTableS1_proc.txt'),
     ('Hs1', '../../data/12-01-start/human Cb660_wan110525_proc.txt'),
     ('Hs2', '../../data/12-01-start/human G166_wan110519_proc.txt'),
     ('Hs3', '../../data/12-01-start/human G166_wan110402_fix_proc.txt'),
     ('Mm0', '../../data/12-01-start/mES_wan110630_proc.txt'),
     ('Dm0', '../../data/12-01-start/Drosophila_embryo_uniprotfasta_wan110728_proc.txt'),
     ('Dm1', '../../data/12-01-start/Drosophila_SL2_uniprotfasta_wan110722_proc.txt'),
     ('Ce0', '../../data/12-01-start/C elegans__uniprotfasta_wan110701_proc.txt'),
     ('Ce1', '../../data/12-01-start/C elegans_uniprotfasta_wan110427_proc.txt'),
     ('Sp0', '../../data/12-01-start/sea urchin hatchedblastula_wan110819_proc.txt'),
     ('Sp1', '../../data/12-01-start/sea urchin unfertlized_wan110902_proc.txt'),
     ('Sp2', '../../data/12-01-start/sea urchin 5minpostunfertlized_wan110908_proc.txt'),
     ('Sp3', '../../data/12-01-start/sea urchin unfertlized_wan110330_proc.txt'),
     ('Sp4', '../../data/12-01-start/sea urchin 2cellculturecleaved_wan110827_proc.txt'),
     ('Dd', '../../data/12-01-start/Dicty discoideum_wan111011_proc.txt')
     ]

def plot_corrs(multi_elut):
    for n in multi_elut.names:
        pl.cluster(multi_elut.eluts[n].corr)

def calc_corrs():
    multi_elut = le.load_multi(elut_files)
    for n in multi_elut.names:
        print 'correlating', n
        elut = multi_elut.multi_elut[n]
        elut.corr = ml.coelution(elut)
    return multi_elut
