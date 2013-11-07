"""
hyppy.py: Python implementation of WFG Hypervolume algorithm
Copyright (C) 2013  Matthew Woodruff, Jon Herman

Original C WFG implementation (C) 2012 Lucas Bradstreet and others

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

----
This implementation assumes minimization.
"""

"""
wfg(pl):
    return sum {exclhv(pl, k) | k in {1 .. |pl|}}

exclhv(pl, k):
    return inclhv(pl[k]) - wfg(nds(limitset(pl, k)))

inclhv(p):
    return product {|p[j] - refPoint[j]| | j in {1 .. n}}

limitset(pl, k):
    for i = 1 to |pl| - k
        for j = 1 to n
            ql[i][j] = worse(pl[k][j], pl[k+i][j])
    return ql

nds(pl) returns the non-dominated subset of pl

"""
import argparse
import sys

class WFG(object):
    def __init__(self, refpoint):
        """
        create a hypervolume-computing object with the given 
        reference point
        """
        self.refpoint = refpoint

    def wfg(self, front):
        """
        return the hypervolume of a front
        front: a list of points

        wfg(pl):
            return sum {exclhv(pl, k) | k in {1 .. |pl|}}
        """
        return sum(self.exclusive(front, index)
                   for index in range(len(front)))

    def exclusive(self, front, index):
        """
        return the exclusive hypervolume of a point relative to a front
        front: a list of points
        index: index of a point in the front

        exclhv(pl, k):
            return inclhv(pl[k]) - wfg(nds(limitset(pl, k)))
        """
        limitset=verboselimitset
        return self.inclusive(front[index])\
               - self.wfg(nds(limitset(front, index)))

    def inclusive(self, point):
        """
        return the product of the difference between all objectives

        inclhv(p):
            return product {|p[j] - refPoint[j]| | j in {1 .. n}}
        and a reference point
        """
        offset = [p-r for (p,r) in zip(point, self.refpoint)]
        volume = 1
        for val in offset:
            volume *= val
        return abs(volume)

def verboselimitset(front, index):
    result = limitset(front, index)
    print "front:"
    for ii in range(len(front)):
        if index == ii:
            print "*", front[ii]
        else:
            print " ", front[ii]
    print "limit set:"
    for row in result:
        print row
    return result

def limitset(front, index):
    """
    return the remainder of the front as limited by the point
    front: a list of points
    index: an index into the list of points

    limitset(pl, k):
        for i = 1 to |pl| - k
            for j = 1 to n
                ql[i][j] = worse(pl[k][j], pl[k+i][j])
        return ql

    """
    return [[max(p,q) for (p,q) in zip(front[j], front[j+index+1])]
            for j in range(len(front)-index-1)]

def nds(front):
    """
    return the nondominated solutions from a set of points
    """
    print "sorting front of size {0}".format(len(front))
    if len(front) <= 1:
        return front
    print "really sorting front of size {0}".format(len(front))
    archive = []

    for row in front:
        asize = len(archive)
        ai = -1
        while ai < asize - 1:
            ai += 1
            adominate = False
            sdominate = False
            nondominate = False
            for arc, sol in zip(archive[ai], row):
                if arc < sol:
                    adominate = True
                    if sdominate:
                        nondominate = True
                        break # stop comparing objectives
                elif arc > sol:
                    sdominate = True
                    if adominate:
                        nondominate = True
                        break # stop comparing objectives
            if nondominate:
                continue # compare next archive solution
            if adominate:
                break    # row is dominated
            if sdominate:
                archive.pop(ai)
                ai -= 1
                asize -= 1
                continue # compare next archive solution
        # if the solution made it all the way through, keep it
        archive.append(row)
    return archive

def linesof(fp):
    """
    yield the lines read from a file
    """
    line = fp.readline()
    while line != "":
        yield line
        line = fp.readline()

def rowsof(lines, sep=" "):
    """
    take a stream of lines and yield a stream of rows
    """
    for line in lines:
        yield [float(x) for x in line.split(sep)]

def onlydata(lines):
    """
    yield lines that are data, break if not
    """
    line = None
    for line in lines:
        if not line.startswith("#"):
            break

    if line is not None:
        yield line

    for line in lines:
        if line.startswith("#"):
            break
        else:
            yield line

def tables_in(lines):
    """
    yield tables from lines
    """
    while True:
        table = list(rowsof(onlydata(lines)))
        if len(table) > 0:
            yield table
        else:
            break


def cli(argv):
    tables = tables_in(linesof(sys.stdin))
    for table in tables:
        print table
        wfg = WFG([0] * len(table[0]))
        try:
            hv = wfg.wfg(table)
            print hv
        except TooDeep:
            pass

if __name__ == "__main__":
    cli(sys.argv)
