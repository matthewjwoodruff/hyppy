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

import argparse
import sys

class HyppyError(Exception): pass

class WFG(object):
    def __init__(self, refpoint):
        """
        create a hypervolume-computing object with the given 
        reference point
        """
        self.refpoint = refpoint
#        self.exclusive = self.verbose_exclusive

    def iterative(self, table):
        stack = [] # tuples: front, lfront, index, inclusive HV, list of exclusive HVs
        front = table
        lfront = len(table)
        index = 0
        excl= []
        depth = 1
        hv = 0

        while depth > 0:
            if index >= lfront: # all points have been processed
                hv = sum(excl)
                depth -= 1
                if depth > 0:
                    front, lfront, index, incl, excl = stack.pop()
                    excl.append(incl - hv)
                    index += 1
            else:
                point = front[index]
                incl = self.inclusive(point)
                limset = limitset(front, index)
                lls = len(limset)
                if lls == 0:
                    excl.append(incl)
                    index += 1
                else:
                    if lls > 1:
                        limset = nds(limset)
                        lls = len(limset)
                    stack.append((front, lfront, index, incl, excl))
                    front = limset
                    lfront = lls
                    index = 0
                    excl = []
                    depth += 1
        return hv

    def wfg(self, front):
        """
        return the hypervolume of a front
        front: a list of points

        wfg(pl):
            return sum {exclhv(pl, k) | k in {1 .. |pl|}}
        """
        return sum(self.exclusive(front, index)
                   for index in range(len(front)))

    def verbose_exclusive(self, front, index):
        """ verbose version of exclusive for debugging """
        incl = self.inclusive(front[index])
        ls = limitset(front, index)
        ls_hv = self.wfg(nds(ls))
        excl = incl - ls_hv

        print("\nfront {0}".format(len(front)))
        for ii in range(len(front)):
            if index == ii:
                print("* {0}".format(front[ii]))
            else:
                print("  {0}".format(front[ii]))
        print("limit set {0}".format(len(ls)))
        for row in ls:
            print("  {0}".format(row))
        print("inclusive hypervolume of front[{0}]: {1}".format(index, incl))
        print("hypervolume of limit set: {0}".format(ls_hv))
        print("exclusive hypervolume of front[{0}]: {1}".format(index, excl))
        return excl

    def exclusive(self, front, index):
        """
        return the exclusive hypervolume of a point relative to a front
        front: a list of points
        index: index of a point in the front

        exclhv(pl, k):
            return inclhv(pl[k]) - wfg(nds(limitset(pl, k)))
        """
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
    """
    for debugging: print out the limit set.
    Use this by setting limitset=verboselimitset in the 
    function where you call limitset.  This could break
    everything.  Use at your own risk.
    """
    result = limitset(front, index)
    print("front:")
    for ii in range(len(front)):
        if index == ii:
            print("* {0}".format(front[ii]))
        else:
            print("  {0}".format(front[ii]))
    print("limit set:")
    for row in result:
        print(row)
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
    return [[max(p,q) for (p,q) in zip(front[index], front[j+index+1])]
            for j in range(len(front)-index-1)]

def nds(front):
    """
    return the nondominated solutions from a set of points
    """
    fsize = len(front)
    if fsize < 2:
        return front
    nobj = len(front[0])
    archive = [False] * fsize

    for sol in xrange(fsize):
        for arc in xrange(fsize):
            if not archive[arc]:
                continue
            adominate = False
            sdominate = False
            nondominate = False
            for ii in xrange(nobj):
                if front[arc][ii] < front[sol][ii]:
                    adominate = True
                    if sdominate:
                        nondominate = True
                        break # stop comparing objectives
                elif front[arc][ii] > front[sol][ii]:
                    sdominate = True
                    if adominate:
                        nondominate = True
                        break # stop comparing objectives
            if nondominate:
                continue # compare next archive solution
            if adominate:
                break    # row is dominated
            if sdominate:
                archive[arc] = False
                continue # compare next archive solution
        # if the solution made it all the way through, keep it
        archive[sol] = True
    return [front[sol] for sol in xrange(fsize) if archive[sol]]

def verbose_hv_of(tables):
    """ yield hypervolume of each table """
    for table in tables:
        print("table {0}".format(len(table)))
        for row in table:
            print(row)
        referencepoint = [0] * len(table[0])
        wfg = WFG(referencepoint)
        yield wfg.wfg(table)

def hv_of(tables):
    """ yield hypervolume of each table """
    for table in tables:
#        table.sort(reverse=True) # in place
        referencepoint = [0] * len(table[0])
        wfg = WFG(referencepoint)
#        yield wfg.wfg(table)
        yield wfg.iterative(table)

def linesof(fp):
    """
    yield the lines read from a file, with line number
    """
    line = fp.readline()
    counter = 0
    while line != "":
        counter += 1
        yield (counter, line)
        line = fp.readline()

def rowsof(lines, delim=" "):
    """
    take a stream of lines and yield a stream of rows
    """
    number = None
    for number, line in lines:
        yield (number, [x for x in line.split(delim)])

def onlydata(lines, ts):
    """
    yield lines that are data, break if not
    """
    line = None
    for number, line in lines:
        if not line.startswith(ts):
            break

    if line is not None:
        yield (number, line)

    for number, line in lines:
        if line.startswith(ts):
            break
        else:
            yield (number, line)

def objectives_in(table, objectives):
    """ yield the objectives in each row """
    number = None
    try:
        if objectives is not None:
            for (number, row) in table:
                obj = [float(row[i]) for i in objectives]
                yield (number, obj)
        else:
            for (number, row) in table:
                obj = [float(x) for x in row]
                yield (number, obj)
    except IndexError as ie:
        msg = "Unable to index objectives on line {0}, error: {1}"
        raise HyppyError(msg.format(number, ie))
    except ValueError as ve:
        msg = "Unable to convert value to float on line {0}, error: {1}"
        raise HyppyError(msg.format(number, ve))
    except TypeError as te:
        msg = "Tried to index with {0} on line {1} of input, error: {2}"
        raise HyppyError(msg.format(objectives, number, te))

def maximize(table, indices=None):
    """
    yield rows of the table with specified indices negated
    table: rows of objectives
    indices: indices into the table. None means maximize all
    """
    number = None
    if indices is None:
        for (number, row) in table:
            yield (number, [-1.0 * x for x in row])
    else:
        try:
            for (number, row) in table:
                newrow = []
                for ii in range(len(row)):
                    if ii in indices:
                        newrow.append(-1.0 * row[ii])
                    else:
                        newrow.append(row[ii])
                yield (number, newrow)
        except IndexError as ie:
            msg = "could not find all maximization objectives on line {0}, error: {1}"
            raise HyppyError(msg.format(number, ie))

def tables_in(lines, **kwargs):
    """ yield tables from lines """
    sep = kwargs.get("separator", None)
    objectives = kwargs.get("objectives", None)
    delim = kwargs.get("delimiter", " ")
    maximize_all = kwargs.get("maximize_all", False)
    max_indices = kwargs.get("maximize", None)

    while True:
        if sep is not None:
            data = onlydata(lines, sep)
        else:
            data = lines
        table = rowsof(data, delim)
        obj = objectives_in(table, objectives)

        if maximize_all:
            obj = maximize(obj)
        elif max_indices is not None:
            obj = maximize(obj, max_indices)

        table = [oo for number, oo in obj]
        if len(table) > 0:
            yield table
        else:
            break

def rerange(intranges):
    """ convert a set of intranges into a list of integers """
    if intranges is None:
        return None
    thelist = []
    for therange in intranges:
        thelist.extend(therange)
    return thelist

def intrange(arg):
    """ convert a command-line argument to a list of integers """
    acceptable_chars = [str(x) for x in range(10)]
    acceptable_chars.append("-")

    partial = []
    first = None

    msg = "Could not convert {0} to index range.".format(arg)
    err = TypeError(msg)

    for char in arg:
        if char not in acceptable_chars:
            raise err
        if char == "-":
            if len(partial) == 0:
                raise err
            elif first is None:
                first = int("".join(partial))
                partial = []
            else: # this means there's a second -, which is not ok
                raise err
        else:
            partial.append(char)

    second = None
    if first is None:
        first = int("".join(partial))
    elif len(partial) == 0:
        raise err
    else:
        second = int("".join(partial))

    if second is None:
        return [first]
    elif second - first >= 0:
        return range(first, second+1)
    else:
        return range(first, second-1, -1)

def argparser(name):
    """ create an argument parser """
    parser = argparse.ArgumentParser(name)
    parser.add_argument("inputfile", type=argparse.FileType("r"),
                        help="file containing sets for which "\
                             "to compute hypervolume")
    parser.add_argument('-o', '--objectives', type=intrange, nargs='+',
                        help='objective columns (zero-indexed)')
    parser.add_argument('-m', '--maximize', type=intrange, nargs='+',
                        help='objective columns to maximize')
    parser.add_argument('-M', '--maximize-all', action="store_true",
                        help='maximize all objectives')
    parser.add_argument("--reverse-column-indices", action='store_true',
                        default=False, help='Reverse the order of column '\
                        'indices.  May be useful if your objectives are '\
                        'at the end of a row of unknown length.  Make sure '\
                        '-e and -m are consistent with the order you '\
                        'specify.')

    delimiters = parser.add_mutually_exclusive_group()
    delimiters.add_argument('-d', '--delimiter', type=str, default=' ',
                        help='input column delimiter, default to space (" ")')
    delimiters.add_argument('--tabs', action="store_true",
                        help="use tabs as delimiter")

    parser.add_argument("-s", "--table-separator",
                        help="character used to separate tables in the input "\
                             "file, default is '#'",
                        default = "#")
    return parser.parse_args

def postprocess(parse):
    """
    return a function that parses argv and does postprocessing
    parse: a function that parses argv
    """
    def postprocd(argv):
        """ the function that parses argv and does postprocessing """
        args = parse(argv)
        args.objectives = rerange(args.objectives)
        args.maximize = rerange(args.maximize)
        if args.tabs:
            args.delimiter = "\t"
        if args.reverse_column_indices:
            if args.objectives is not None:
                args.objectives = [-1 - ob for ob in args.objectives]
            if args.maximize is not None:
                args.maximize = [-1 - ob for ob in args.maximize]
        if args.maximize is not None and args.objectives is not None:
            try:
                # transform to an index into the objectives
                args.maximize = [args.objectives.index(i) for i in args.maximize]
            except ValueError:
                raise HyppyError("cannot maximize an objective that does not exist")
        return args
    return postprocd

def cli(argv):
    """ command-line interface to hyppy """
#    hv_of = verbose_hv_of
    parse = postprocess(argparser(argv.pop(0)))
    args = parse(argv)

    lines = linesof(args.inputfile)
    tables = tables_in(lines,
                       objectives = args.objectives,
                       maximize = args.maximize,
                       maximize_all = args.maximize_all,
                       delimiter = args.delimiter,
                       separator = args.table_separator)

    for hv in hv_of(tables):
        print(hv)

if __name__ == "__main__":
    cli(sys.argv)
