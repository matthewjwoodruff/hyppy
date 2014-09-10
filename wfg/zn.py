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
    nadir = [0.0]*nobj
    zenith = [0.0]*nobj
    vol = 0.0

    for row in rows:
        for i, val in enumerate(row):
            nadir[i] = min([nadir[i], val])
            zenith[i] = max([zenith[i], val])

    add_volume = 1.0
    for z in zenith:
        add_volume *= z
    sub_volume = 1.0
    for z, n in zip(zenith, nadir):
        sub_volume *= (z-n)
    vol += add_volume
    vol -= sub_volume

    return vol
