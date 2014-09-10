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
    if len(rows) == 0:
        return 0.0
    nobj = len(rows[0])
    nadir = [float("inf")]*nobj
    zenith = [-float("inf")]*nobj
    vol = 0.0

    for row in rows:
        for i, val in enumerate(row):
            nadir[i] = min([nadir[i], val])
            zenith[i] = max([zenith[i], val])
    print("nadir {0}".format(nadir))
    print("zenith {0}".format(zenith))

    add_volume = 1.0
    for z in zenith:
        add_volume *= z
    print("add {0}".format(add_volume))
    sub_volume = 1.0
    for z, n in zip(zenith, nadir):
        sub_volume *= (z-n)
    print("sub {0}".format(sub_volume))
    vol += add_volume
    vol -= sub_volume

    # translate nadir to origin
    transformed = []
    for row in rows:
        transformed.append([r-n for r, n in zip(row, nadir)])

    # discard extreme points
    transformed = [
            r for r in transformed 
            if not any([x <= 0 for x in r])]

    # compute hv on the remaining points
    vol += zn(transformed)

    # this could *easily* be made tail-recursive by passing vol
    return vol
