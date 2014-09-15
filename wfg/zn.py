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
    if len(rows) == 0:
        return 0
    global calls_at_level
    calls_at_level[len(rows[0])] = calls_at_level.get(len(rows[0]), 0)
    calls_at_level[len(rows[0])] += 1

    nobj = len(rows[0])
    if nobj == 0:
        return 1.0
    vol = 0.0

    while len(rows) > 0:
        # find hypervolume one dimension down for each silhouette
        for axis in range(nobj):
            if len(rows) == 0:
                break
            nadir = float("inf")
            nadir_index = None
            secondary_nadir = float("inf")

            for rowindex, row in enumerate(rows):
                if row[axis] < nadir:
                    secondary_nadir = nadir
                    nadir = row[axis]
                    nadir_index = rowindex
                elif row[axis] == nadir:
                    nadir_index = None

            sil = []
            carry_along = []

            offset = nadir
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

            # Now if there's exactly one nadir point and it's in the silhouette, 
            # we can compute its one-down hypervolume easily, do another step, 
            # and remove it from the set.  Except not really.
            # The secondary nadir is not sufficient. We also need the one-down
            # hypervolume of points that are one-down dominated only by the nadir
            # point we're removing.
            if nadir in carry_along and secondary_nadir != float("inf") and False:
                remove_point = sil[carry_along.index(nadir)]
                offset2 = secondary_nadir - nadir
                subtract_box = sum([x**2 for x in remove_point])
                vol += offset2 * down_one
                vol -= offset2 * subtract_box
                #print("subtracted!")
                #exit()
                offset += offset2

            newrows = []
            for row in rows:
                if row[axis] > offset:
                    row[axis] -= offset
                    newrows.append(row)
            rows = newrows

    return vol

def zn(rows):
    global calls_at_level
    calls_at_level = dict()
    hv = _zn(rows)
    print(calls_at_level)
    return hv
