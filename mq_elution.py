import re
import elution as el
import utils as ut

def mq2elut(fname, quant='iBAQ'):
    lines = [l for l in ut.load_tab_file(fname)]
    # want eg 'iBAQ WAN...', not 'iBAQ L WAN...'
    inds = [i for i,val in enumerate(lines[0]) 
            if re.match('^%s\s\w{2}' % quant,val) is not None]
    #prots = [[p.split()[0][1:] for p in ps.split(';')] 
            #for ps in [l[0] for l in lines[1:]]]
    # for now just using the "majority protein"
    prots = [p.split()[0][1:] for p in [l[1] for l in lines[1:]]]
    output = [[lines[0][0]] + [lines[0][i] for i in inds]] + \
            [[p] + [l[i] for i in inds] for p,l in zip(prots, lines[1:])]
    ut.write_tab_file(output, ut.pre_ext(fname, '_mq_%s' % quant))



