from __future__ import division
import utils as ut
from collections import defaultdict
import itertools as it

def collapse_alternatives(ppis):
    """
    Simply go through the whole set of interactions, and whenever there are two
    nodes with identical partners, and scores are all within 20% or 0.01,
    collapse them to a single node.
    """
    nodes = set([p[0] for p in ppis]+[p[1] for p in ppis])
    d_ints = defaultdict(set)
    for node1,node2,score in ppis:
        d_ints[node1].add((node2,score)); d_ints[node2].add((node1,score))
    alt_groups = []
    for i,(node1,ints1) in enumerate(d_ints.items()):
        # first go through the existing collections
        found_group=False
        for group in alt_groups:
            anode1,an1ints = group[0]
            if (len(ints1)==len(an1ints) and
                    set(ut.i0(ints1))==set(ut.i0(an1ints))):
                group.append((node1,ints1))
                found_group=True
                break
        if not found_group:
            alt_groups.append([(node1,ints1)])
    return alt_groups

