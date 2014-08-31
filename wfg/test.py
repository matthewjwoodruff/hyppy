from ctypes import CDLL, c_double

reduced = []
nobj=0
npoints=0
with open('./reduced', 'r') as fp:
    for line in fp:
        if "#" in line:
            continue
        row = [float(x) for x in line.split()]
        if nobj > 0 and len(row) != nobj:
            raise Exception("bad!")
        nobj = len(row)
        npoints += 1
        reduced.extend(row)

print("reduced is {0} long, has {1} objectives and {2} points".format(len(reduced), nobj, npoints))

wfg2 = CDLL('./libwfg2_debug.so')
hv = wfg2.hypervolume
hv.restype = c_double
Dub = c_double * (nobj * npoints)
dub = Dub(*reduced)
print(hv(nobj, npoints, dub))

