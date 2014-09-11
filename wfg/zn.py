"""
zn.py
Fast hypervolume computation.
Maximization relative to origin.
Matthew Woodruff
2014-09-10
LGPL
"""
import copy
current_depth = -1

def zn(rows):
    """
    rows: the set of rows for which to compute hypervolume
    Assume maximization relative to origin.
    """
    global current_depth
    current_depth += 1
    indent = "".join(["  "]*current_depth)
    print("".join([indent, "call at depth {0}, rows {1}".format(current_depth, rows)]))
    nobj = len(rows[0])
    if nobj == 0:
        print("".join([indent,"return 1.0 at depth {0}, rows {1}".format(current_depth, rows)]))
        current_depth -= 1
        return 1.0
    nadir = [float("inf")]*nobj
    zenith = [-float("inf")]*nobj
    zenith_contributors = [nadir] * nobj
    vol = 0.0

    # discard points that are at zero
    rows = [
            r for r in rows 
            if not any([x <= 0 for x in r])]

    if len(rows) == 0:
        print("".join([indent,"return 0.0 at depth {0}, rows {1}".format(current_depth, rows)]))
        current_depth -= 1
        return 0.0

    for row in rows:
        for i, val in enumerate(row):
            nadir[i] = min((nadir[i], val))
            if val > zenith[i]:
                zenith[i] = val
                zenith_contributors[i] = copy.copy(row)
            elif val == zenith[i]:
                if any(r > zc for r, zc in zip(row, zenith_contributors[i])):
                    zenith_contributors[i] = copy.copy(row)

    # find hypervolume one dimension down for each silhouette
    for axis in range(nobj):
        offset = nadir[axis]
        sil = []
        for row in zenith_contributors:
            sil.append([row[i] for i in range(nobj) if i != axis])
        sil = [list(row) for row in list(set([tuple(row) for row in sil]))]
        down_one = zn(sil)
        print("".join(
            [indent,"silhouette for axis {0} is {1}".format(axis, sil)]))
        print("".join(
            [indent, 
             "down one hypervolume is {0} with offset {1} for a contribution of {2}".format(
                           down_one, offset, offset*down_one)]))
        vol += offset * down_one
        print("".join([indent, "for a running volume of {0}".format(vol)]))
        for row in zenith_contributors:
            row[axis] -= offset

    nadir_hypervolume = 1.0
    for x in nadir:
        nadir_hypervolume *= x
    print("".join(
        [indent, 
         "depth {0} nadir HV is {1}".format(current_depth, nadir_hypervolume)]))
    
    # translate nadir to origin
    transformed = []
    for row in rows:
        transformed.append([r-n for r, n in zip(row, nadir)])

    # compute hv on the translated points
    vol += zn(transformed)

    print("".join(
        [indent, 
         "return at depth {0}, rows {1}".format(current_depth, rows)]))
    current_depth -= 1
    return vol
