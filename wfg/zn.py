"""
zn.py
Fast hypervolume computation.
Maximization relative to origin.
Matthew Woodruff
2014-09-10
LGPL
"""

def zn(rows):
    """
    rows: the set of rows for which to compute hypervolume
    Assume maximization relative to origin.
    """
    nobj = len(rows[0])
    nadir = [float("inf")]*nobj
    zenith = [-float("inf")]*nobj
    vol = 0.0

    # discard points that are at zero
    rows = [
            r for r in rows 
            if not any([x <= 0 for x in r])]

    if len(rows) == 0:
        return 0.0

    for row in rows:
        for i, val in enumerate(row):
            nadir[i] = min([nadir[i], val])
            zenith[i] = max([zenith[i], val])

    # add nadir volume
    nvol = 1.0
    for n in nadir:
        nvol *= n
    vol += nvol

    # add extra chunks on the sides
    for i, z in enumerate(zenith):
        point = [n for n in nadir]
        point[i] = z - nadir[i]
        chunk = 1.0
        for x in point:
            chunk *= x
        vol += chunk

    # translate nadir to origin
    transformed = []
    for row in rows:
        transformed.append([r-n for r, n in zip(row, nadir)])

    # compute hv on the remaining points
    vol += zn(transformed)

    # this could *easily* be made tail-recursive by passing vol
    return vol
