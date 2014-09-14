"""
zn.py
Fast hypervolume computation.
Maximization relative to origin.
Matthew Woodruff
2014-09-10
LGPL
"""
import copy

calls_at_level = dict()
step = 0

def sort(rows):
    """
    true Pareto sort assuming maximization
    """
    archive = []
    asize = 0

    for row in rows:
        ai = -1
        adominate = False # archive dominates
        sdominate = False # solution dominates
        nondominate = False # neither dominates
        while ai < asize - 1:
            ai += 1
            adominate = False
            sdominate = False
            nondominate = False
            aobj = archive[ai]
            for ar, ro in zip(aobj, row):
                if ar > ro:
                    adominate = True
                    if sdominate:
                        nondominate = True
                        break # for
                elif ar < ro:
                    sdominate = True
                    if adominate:
                        nondominate = True
                        break # for
            if nondominate:
                continue
            if adominate:
                break
            if sdominate:
                archive.pop(ai)
                ai -= 1
                asize -= 1
                continue 
        if nondominate or sdominate or not adominate:
            archive.append(row)
            asize += 1
    return archive

def _zn(rows):
    """
    rows: the set of rows for which to compute hypervolume
    Assume maximization relative to origin.
    """
    global calls_at_level
    global step
    level = len(rows[0])
    calls_at_level[level] = calls_at_level.get(level, 0)
    calls_at_level[level] += 1

    nobj = len(rows[0])
    if nobj == 0:
        return 1.0
    nadir = [float("inf")]*nobj
    zenith = [-float("inf")]*nobj
    nadir_contributors = [None]*nobj # indices of nadir components
    vol = 0.0

    # discard points that are at zero
    rows = [
            r for r in rows 
            if not any([x <= 0 for x in r])]

    if len(rows) == 0:
        return 0.0

    for row in rows:
        for i, val in enumerate(row):
            nadir[i] = min((nadir[i], val))

    # find hypervolume one dimension down for each silhouette
    for axis in range(nobj):
        offset = nadir[axis]
        if nobj > 1:
            sil = []
            for row in rows:
                sil.append([row[i] for i in range(nobj) if i != axis])
            sil = sort(sil)
            down_one = _zn(sil)
        else:
            down_one = 1.0
        vol += offset * down_one
        for row in rows:
            row[axis] -= offset

    # compute hv on the translated points
    step += 1 # count the steps
    vol += _zn(rows)

    return vol

def zn(rows):
    global calls_at_level
    global step
    hv = _zn(rows)
    print(calls_at_level)
    print("{0} steps".format(step))
    calls_at_level = dict()
    step = 0
    return hv
