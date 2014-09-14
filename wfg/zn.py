"""
zn.py
Fast hypervolume computation.
Maximization relative to origin.
Matthew Woodruff
2014-09-10
LGPL
"""
calls_at_level = dict()

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
    level = len(rows[0])
    calls_at_level[level] = calls_at_level.get(level, 0)
    calls_at_level[level] += 1

    nobj = len(rows[0])
    if nobj == 0:
        return 1.0
    vol = 0.0

    while len(rows) > 0:
        nadir = [float("inf")]*nobj
        #nadir_contributors = [None]*nobj # indices of nadir components
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
                if nobj == 2 and len(sil) == 1:
                    down_one = sil[0][0] # one dimensional hv is easy
                else:
                    down_one = _zn(sil)
            else:
                down_one = 1.0
            vol += offset * down_one
            for row in rows: # translates to nadir
                row[axis] -= offset

        # discard points that are at zero
        rows = [
                r for r in rows 
                if not any([x <= 0 for x in r])]

    return vol

def zn(rows):
    global calls_at_level
    hv = _zn(rows)
    print(calls_at_level)
    calls_at_level = dict()
    return hv
