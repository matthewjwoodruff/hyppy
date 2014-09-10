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
    zenith_contributors = [nadir] * nobj
    vol = 0.0

    # discard points that are at zero
    rows = [
            r for r in rows 
            if not any([x <= 0 for x in r])]

    if len(rows) == 0:
        return 0.0

    for row in rows:
        for i, val in enumerate(row):
            if val > zenith[i]:
                zenith[i] = val
                zenith_contributors[i] = row
            elif val == zenith[i]:
                if any(r > zc for r, zc in zip(row, zenith_contributors[i])):
                    zenith_contributors[i] = row
        print("row {0}, zc {1}".format(row, zenith_contributors))

    for zc in zenith_contributors:
        for i, val in enumerate(zc):
            nadir[i] = min([nadir[i], val])

    for axis in range(nobj-1):
        # step along axis, add volumes from next axis, with other
        # dimensions clamped to nadir
        inorder = sorted(zenith_contributors, key=lambda row: row[axis])
        inorder = [[y for y in x] for x in inorder]
        relative_to = [0.0] * nobj

        for i in range(axis): # axes already visited
            relative_to[i] = nadir[i]
        print("inorder {0}, relative to {1}".format(inorder, relative_to))

        while len(inorder) > 0:
            point = inorder[0]

            # axis+1 gets the best value remaining
            #point[axis+1] = max([x[axis+1] for x in inorder])

            # higher axes get nadir value
            for i in range(axis+2, nobj):
                point[i] = nadir[i]

            # compute volume
            chunk = 1.0
            for x, r in zip(point, relative_to):
                chunk *= (x-r)
            vol += chunk

            print("adding {0} for {1} relative to {2}".format(chunk, point, relative_to))

            # update point relative to
            relative_to[axis] = point[axis]

            # remove point
            inorder.pop(0)

    # translate nadir to origin
    transformed = []
    for row in rows:
        transformed.append([r-n for r, n in zip(row, nadir)])

    # compute hv on the remaining points
    vol += zn(transformed)

    # this could *easily* be made tail-recursive by passing vol
    return vol
