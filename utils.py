from __future__ import division
import cPickle
import sys
import os
#from scipy import array, random # causing errors weirdly on mac air
#import scipy
import random
import operator
import itertools
import string
import warnings
import numpy as np



######################################################################
## PYTHON OBJECT SAVE/LOAD
######################################################################

def savepy(obj,fname,safe=False):
    # If safe is true, we save obj to a temporary file first, then mv the file
    # to its final destination
    if safe:
        save(obj, fname+'.partial', safe=False)
        shutil.move(fname+'.partial', fname)
        return
    # Because reloading a module thwarts pickling (the new class is not the same
    # as the old)
    def maybe_reload(obj):
        if hasattr(obj, '_pickle_reload'):
            obj = obj._pickle_reload()
        return obj
    obj = maybe_reload(obj)
    cPickle.dump(obj, file(fname, 'wb'), protocol=2)

def loadpy(fname):
    obj = cPickle.load(file(fname, 'rb'))
    # if isinstance(obj, bayc.World):
    #     obj.update_version()
    if hasattr(obj, '_pickle_update'):
        # This function is called on loading, and should return the new
        # object, or self if the object has updated itself.
        obj = obj._pickle_update()
    try:
        obj.filename = fname
    except AttributeError:
        pass # Some objects don't have a filename
    return obj


########################################################################
## COLLECTIONS and math functions
########################################################################

#t_array = type(array([1])) # because type(array([1])) != array

def all_same(f, bag):
    v = f(bag[0])
    for x in bag[1:]:
        if f(x) != v: return False
    return True

def bin(list, binsize):
    nbins = int(np.ceil(len(list)/binsize))
    return [list[i*binsize:(i+1)*binsize] for i in range(nbins)]

def column_totals(mat):
    return np.array(mat.sum(axis=0)).flatten()

def column_uniques(mat):
    return np.array((mat > 0).sum(axis=0)).flatten()

def every(pred, bag):
    # Like the Common Lisp EVERY
    for x in bag:
        if not pred(x): return False
    return True

def fremove(el, l): # functional remove
    l2 = l[:]
    l2.remove(el)
    return l2

def fsort(l, *args, **kwargs):
    assert(type(l) == list)
    l2 = l[:]
    l2.sort(*args, **kwargs)
    return l2

def hyperg(k, N, m, n):
    exact=1
    return comb(m,k,exact) * comb(N-m,n-k,exact)/ comb(N,n,exact)

def hyperg_cum2(c, N, m, n):
    # Taking 1-sum_from_0_to_c is much, much faster than summing from c
    # to min(m,n) as is done in the paper
    # FIXME: We lose all of our precision substracting from 1
    return 1-sum([hyperg(k,N,m,n) for k in range(0, c)])

def hyperg_cum(c, N, m, n):
    if c == 0:
        return 1
    else:
        s = sum([hyperg(k,int(N),int(m),int(n)) for k in range(c, min(m,n)+1)])
        # Sanity check. Sadly, it has failed before. Calling
        # hyperg(2,5022,22,17) yields a negative number when the 22 is a
        # numpy.int64, but only if comb() is exact (see help(comb))
        assert(s>0)
        return s

def idpr(x): #useful sometimes for printing nested in expressions
    print x
    return x

def least(f, seq, comp=operator.lt):
    # Returns the x in seq for which f(x) is smallest
    l = None
    lval = None
    for el in seq:
        v = f(el)
        if (lval and comp(v, lval)) or (not lval):
            l = el
            lval = v
    return l

def most(f, seq): return least(lambda x: -f(x), seq, comp=operator.gt)

def only(l):
    assert(len(l)==1)
    return l[0]

def onlys(s): # for set
    return only(tuple(s))

def printnow(s):
    print s
    sys.stdout.flush()

def ravg(l):
    return rsum(l) / len(l)

def rsum(l): #reduce sum
    # This is exactly like the Python built-in sum, but scipy overloads it
    # with its own sum function.
    return reduce(operator.add, l)

# def rnd(a,b=None):
#     # excludes b. we return an int if a is an int, a float if it's a float
#     if isinstance(a,int):
#         if b: assert(isinstance(b,int))
#         if b is None:
#             b = a
#             a = 0
#         return random.randint(a, b-1)
#     else:
#         if b is None:
#             b = a
#             a = 0
#         return random.uniform(a,b)

def sample_wr(pop, k):
    n = len(pop)
    _random, _int = random.random, int
    return [pop[_int(_random() * n)] for i in itertools.repeat(None,k)]

def zip_exact(*seqs):
    # Like zip, but generates an error if the seqs are not all the same
    assert(all_same(len, seqs))
    return zip(*seqs)



#####################################################################
##  LOADING DATA
#####################################################################

def pre_ext(fname, preext):
    path, name = os.path.split(fname)
    basename,ext = os.path.splitext(name)
    return os.path.join(path, basename+preext+ext)

def shortname(fname):
    return os.path.splitext(os.path.basename(fname))[0]

def load_dict_flat(fname):
    # insted use dict(load_tab_file(fname))
    assert 0==1
    pass

def load_tab_file(fname, sep='\t'):
    """ Returns a generator of a list of list
    (assumes each element in a line is tab-delimited)
    """
    def clean(x):
        return x.strip() # remove the newlines
    return (tuple(clean(l).split(sep)) for l in file(fname, 'r'))

def load_dict_sets(fname, lines='pairs'):
    # Returns a key->set(values) mapping
    if lines=='pairs':
        return dict_sets_from_tuples(load_list_of_lists(fname))
    elif lines=='lists':
        return dict([(l[0],set(l[1:])) for l in load_list_of_lists(fname)])

def load_list(fname):
    return [row[0] for row in load_tab_file(fname)]

def load_array(fname):
    return np.loadtxt(fname)

def load_list_of_lists(fname, **kwargs):
    return [l for l in load_tab_file(fname, **kwargs)]

def write_dict_sets_lines(d,fname):
    write_tab_file([[i[0]]+list(i[1]) for i in d.items()],fname)
    
def write_dict_sets_tab(d,fname):
    mylist = []
    for k in d:
        mylist.extend([(k,v) for v in d[k]])
    write_tab_file(mylist,fname)

def write_tab_file(ll, fname, formatter='{0}'):
    # '{0:10g}' for 10 digit precision
    f = file(fname, 'w')
    for l in ll:
        if isinstance(l, type('teststring')):
            f.write(formatter.format(l))
        else:
            for i,x in enumerate(l):
                if i>0: f.write('\t')
                f.write(formatter.format(x))
        f.write('\n')

def print_lol(lol):
    for l in lol:
        print '\t'.join([str(i) for i in flatten(l)])

def print_list(l):
    print '\t'.join([str(i) for i in l])

def flatten(lst):
   out = []
   for sublist in lst:
       if not isinstance(sublist,list):
           out.append(sublist)
       else:
           out.extend(sublist)
   return out



#######################################
# DICT OPERATIONS
#######################################


def dict_remove_if_value(d,remove_val):
    emptys = []
    for key in d:
        if d[key] == remove_val:
            emptys.append(key)
    for k in emptys:
        del(d[k])
    return d

def compose_dict_sets(d1,d2,compose_op=set.union,default=set()):
    # d1 must be key1->set(keys2)
    # and d2 be key2->V
    # compose_op must accept two V
    # default is for d2, it should return a V
    return dict(((k,reduce(compose_op, [d2.get(k2,default) for k2 in v],
                           default)) for (k,v) in d1.items()))

def dict_combine_sets(d1,d2):
    # Assume that the values are sets, and combine them pairwise (union)
    return dict([(k, d1.get(k,set()).union(d2.get(k,set())))
                 for k in list(set(d1.keys()+d2.keys()))])

def dict_inverse_sets(d):
    # If d is a K->set(V) mapping, return a V->set(K) dict
    dout = {}
    for key, values in d.items():
        for v in values:
            dout.setdefault(v,set([])).add(key)
    return dout

def dict_inverse(d):
    dout = {}
    for k, v in d.items():
        dout[v] = k
    return dout

def list_inv_to_dict(lst):
    d = {}
    for index,item in enumerate(lst):
        d[item]=index
    return d

def dict_sets_from_tuples(lot):
    d = {}
    for k,v in [tup for tup in lot if len(tup)>1]:
        d.setdefault(k,set([])).add(v)
    return d

def dict_count_values(d):
    return sum([len(d[k]) for k in d])

def dict_dedup(d):
    # in a dict of sets, for every k:v[i] found, remove v[i]:k
    d_dedup = d.copy()
    for k in d:
        for vi in d[k]:
            if vi in d_dedup:
                if k in d_dedup[vi]:
                    d_dedup[vi].remove(k)
    return d_dedup


##########################################
# Other convenience
##########################################

dir_project = "~/Dropbox/complex"
base_path = os.path.expanduser(dir_project)

def proj_path(pathkey,basename=''):
    if basename:
        return os.path.join(base_path, config()[pathkey], basename)
    else:
        return os.path.join(base_path, config()[pathkey])

def config():
    def str_to_bool(string):
        return (True if string.lower() == 'true' else False if string.lower() ==
            'false' else string)
    conf_path = os.path.join(base_path, 'config.txt')
    dconf = dict([(l.split()[0],str_to_bool(l.split()[1])) for l in load_list(conf_path)])
    return dconf
    
