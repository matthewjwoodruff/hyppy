"""
zn.py
Fast hypervolume computation.
Maximization relative to origin.
Matthew Woodruff
2014-09-10
LGPL
"""
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

    nobj = len(rows[0])
    if nobj == 0:
        return 1.0
    vol = 0.0

    while len(rows) > 0:
        nadir = [float("inf")]*nobj
        second_nadir = [float("inf")]*nobj
        nadir_contributors = [None]*nobj # indices of nadir components
        nadir_backup = [None]*nobj # indices of backup nadir components

        if len(rows) == 1:
            vol += sum([rows[0][i]**2 for i in range(nobj)])
            break

        for rowindex, row in enumerate(rows):
            for i, val in enumerate(row):
                if val < nadir[i]:
                    nadir[i] = val
                    nadir_backup[i] = nadir_contributors[i]
                    nadir_contributors[i] = rowindex
                elif val < second_nadir[i]:
                    nadir_backup[i] = rowindex
                elif val == nadir[i]: # there can be only one
                    nadir_contributors[i] = None

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

            # now if there's exactly one nadir point, we can compute its
            # one-down hypervolume easily, do another step, and remove it 
            # from the set
            if nadir_contributors[axis] is not None and nadir_backup[axis] is not None:
                contributor = nadir_contributors[axis]
                backup = nadir_backup[axis]
                offset2 = rows[backup][axis] - rows[contributor][axis]
                vol += offset2 * down_one
                subtract_box = sum([rows[contributor][i]**2 for i in range(nobj) if i != axis])
                vol -= offset2 * subtract_box
                for a in range(nobj): # can't use it twice
                    if nadir_contributors[a] == contributor:
                        nadir_contributors[a] = None
                        nadir_backup[a] = None
                offset += offset2

            for row in rows: # translates to nadir
                row[axis] -= offset

        # discard points that are at zero or less
        rows = [
                r for r in rows 
                if not any([x <= 0 for x in r])]

    return vol

def zn(rows):
    hv = _zn(rows)
    return hv
