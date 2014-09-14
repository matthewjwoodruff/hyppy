"""
zn.py
Fast hypervolume computation.
Maximization relative to origin.
Matthew Woodruff
2014-09-10
LGPL
"""
calls_at_level = dict()

def pareto(rows, carry_along):
    """
    true Pareto sort assuming maximization
    carrys along additional info
    """
    archive = []
    carry_along_archive = []
    asize = 0

    for row, calong in zip(rows, carry_along):
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
                carry_along_archive.pop(ai)
                ai -= 1
                asize -= 1
                continue 
        if nondominate or sdominate or not adominate:
            archive.append(row)
            carry_along_archive.append(calong)
            asize += 1
    return archive, carry_along_archive

def _zn(rows):
    """
    rows: the set of rows for which to compute hypervolume
    Assume maximization relative to origin.
    """

    global calls_at_level
    calls_at_level[len(rows[0])] = calls_at_level.get(len(rows[0]), 0)
    calls_at_level[len(rows[0])] += 1

    nobj = len(rows[0])
    if nobj == 0:
        return 1.0
    vol = 0.0

    while len(rows) > 0:
        nadir = [float("inf")]*nobj

        for rowindex, row in enumerate(rows):
            for i, val in enumerate(row):
                if val < nadir[i]:
                    nadir[i] = val

        # find hypervolume one dimension down for each silhouette
        for axis in range(nobj):
            sil = []
            carry_along = []

            offset = nadir[axis]
            if nobj > 1:
                for row in rows:
                    sil.append([row[i] for i in range(nobj) if i != axis])
                    carry_along.append(row[axis])
                sil, carry_along = pareto(sil, carry_along)
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
            if len(carry_along) > 1:
                least = float("inf")
                second_least = float("inf")
                lindex  = None
                slindex = None
                for idx, calong in enumerate(carry_along):
                    if calong < least:
                        second_least = least
                        slindex = lindex
                        least = calong
                        lindex = idx
                    if calong == least:
                        lindex = None
                        slindex = None
                    elif calong < second_least:
                        second_least = calong
                        slindex = idx
                if lindex is not None and slindex is not None:
                    offset2 = carry_along[slindex] - carry_along[lindex]
                    subtract_box = sum([sil[lindex][i]**2 for i in range(nobj-1)])
                    vol += offset2 * down_one
                    vol -= offset2 * subtract_box
                    offset += offset2

            for row in rows: # translates to nadir
                row[axis] -= offset

        # discard points that are at zero or less
        rows = [
                r for r in rows 
                if not any([x <= 0 for x in r])]

    return vol

def zn(rows):
    global calls_at_level
    calls_at_level = dict()
    hv = _zn(rows)
    print(calls_at_level)
    return hv
