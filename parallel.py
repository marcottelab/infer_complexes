from IPython.parallel import Client, Reference
import os
import operator
import utils as ut

def pcollect(f, items, freduce=operator.add, view=None, njobs=6):
    """
    - f: a function that takes a list of items and returns something to be
      collected.
    - items: the full list of items.
    - freduce: how to collect the items.
    """
    



