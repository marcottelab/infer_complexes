from __future__ import division


def truth_table(positives, negatives, cutoff=0.05):
    """
    Compute a truth table for the scores being over threshold for two sets of
    scores.
    """
    def npassing(scores, cutoff):
        return len([x for x in scores if x >= cutoff])
    TP = npassing(positives, cutoff)
    FP = len(positives) - TP
    FN = npassing(negatives, cutoff)
    TN = len(negatives) - FN
    print "TP FP ; FN TN"
    print TP, FP, ';', FN, TN

