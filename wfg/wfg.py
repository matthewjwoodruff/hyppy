"""
Copyright (C) 2014 Matthew Woodruff

This script is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This script is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this script. If not, see <http://www.gnu.org/licenses/>.
===========================================================
wfg.py
"""

from ctypes import CDLL, c_double
import os

class WFGError(Exception): pass

wfg2 = None
hv = None

def set_wfg2(path='.'):
    global wfg2
    global hv
    wfg2 = CDLL(os.path.join(path, 'libwfg2.so'))
    hv = wfg2.hypervolume
    hv.restype = c_double

def wfg(rows):
    """
    compute hypervolume on a set of minimization objective rows using wfg
    rows (list of lists)
    """
    global hv
    if hv is None:
        set_wfg2()
    nobj = None
    values = []
    for row in rows:
        if nobj is None:
            nobj = len(row)
        if nobj != len(row):
            raise WFGError("rows are uneven lengths")
        values.extend(row)
    Arr = c_double * (nobj * len(rows))
    arr = Arr(*values)
    if(len(values) > 1):
        vol = hv(nobj, len(rows), arr)
        return vol
    else:
        raise WFGError('zero values!!!!!')
