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
metricsystem.py
"""

import argparse
import ctypes
import sys
import math
import pareto
import wfg
import zn

class InputError(Exception): pass
class ReferencePointError(Exception): pass

def rerange(intranges):
    """ 
    Convert a set of intranges into a list of integers.
    Copied from the implementation in pareto.py.
    """
    if intranges is None:
        return None
    thelist = []
    for therange in intranges:
        thelist.extend(therange)
    return thelist

def intrange(arg):
    """ 
    Convert a command-line argument to a list of integers.
    Copied from the implementation in pareto.py.
    """
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

def regex(arg):
    """
    converts arg into a compiled regular expression
    arg (str) a command-line argument
    """
    import re
    # allow the re.error to propagate up --- it's pretty informative
    return re.compile(arg)

def get_args(argv):
    """ 
    Get command-line arguments.  Inspired by and partly copied from the
    implementation in pareto.py.
    """
    prog = argv.pop(0)
    parser = argparse.ArgumentParser(
        prog=prog, description='Metrics Calculation for Pareto Approximations')

    # inputs and output
    parser.add_argument('inputs', type=argparse.FileType('r'), nargs='+',
        help='input filenames, use - for standard input')
    parser.add_argument('--output', type=argparse.FileType('w'),
        default=sys.stdout,
        help='output filename, defaulting to stdout')

    # choosing columns
    parser.add_argument('-o', '--objectives', type=intrange, nargs='+',
        help='objective columns (zero-indexed), default is all')
    parser.add_argument('--objective-column-names', type=str, nargs='+',
        help='Names of objective columns')
    parser.add_argument('-m', '--maximize', type=intrange, nargs='+',
        help='Columns to maximize, default is none. '\
             'All columns for which maximization is specified '\
             'are considered to be objectives. '\
             'Maximized objectives are multiplied by -1.0 and '\
             'treated as minimization objectives.'\
             'The final hypervolume is also corrected by negation '\
             'if there is an '\
             'odd number of maximization objectives.')
    parser.add_argument('-M', '--maximize-all', action='store_true',
        help='maximize all objectives')
    parser.add_argument('--maximize-column-names', type=str, nargs='+',
        help='Names of columns to maximize.  These are automatically '\
             'considered to be objectives')
    parser.add_argument("--reverse-column-indices", action='store_true',
        default=False, help='Reverse the order of column '\
        'indices.  May be useful if your objectives are '\
        'at the end of a row of unknown length.  Make sure '\
        '-e and -m are consistent with the order you '\
        'specify.')

    # epsilons
    parser.add_argument('-e', '--epsilons', type=float, nargs='+',
        help='Epsilons, one per objective.  Specifying epsilons '\
             'causes data to be passed through an epsilon-nondomination '\
             'sort before hypervolume is computed.')
    parser.add_argument('--integer', action='store_true',
        help='Use an integer version of the hypervolume computation, '\
             "based on epsilon box index relative to reference point's "\
             'epsilon box index.  Epsilon box indexes are always '\
             'referenced to 0 so that box boundaries are consistent '\
             'with those used in the sort.  Specifying ths option '\
             'will significantly change the computed '\
             'hypervolumes.  This option has no effect if epsilons '\
             'are not specified.')

    # reference point
    parser.add_argument('-R', '--reference', type=float, nargs='+',
        help='Reference point.  If one value is specified, it is '\
             'duplicated for all objectives.  If multiple values '\
             'are specified, only solution sets with the same '\
             'number of objectives will be evaluated.  '\
             'Reference point must be a nadir.  Any inferior '\
             'points will be clamped to the reference.')
    parser.add_argument('--auto-reference', 
        choices=['zero', 'first', 'full'],
        help='Automatic reference point behavior. '\
             '"zero" means that the origin is used.  '\
             '"first" means that the nadir point for the first '\
             'solution set in the input is used as a reference '\
             'for every solution set. '\
             '"full" means that the nadir point for all sets '\
             'combined is used to compute hypervolume for each. '\
             'The full option means that output will be written twice, '\
             'so that an aborted run still produces useable output. '\
             'The first-pass output will include the nadir point used '\
             'for each set.'
             'Combining auto reference with a specified reference point '\
             'means that the reference point will be the nadir of the '\
             'specified and automatic points.'\
             'If neither this option nor a reference point is '\
             'specified, the "zero" behavior is applied.')

    # scaling
    parser.add_argument('--scale',
        choices=['none', 'epsilon', 'reference'],
        help='Scaling mode.  "none" is the default.  "epsilon" '\
             'is synonymous with "none" unless epsilons are '\
             'specified, in which case all objective values '\
             'are scaled by their epsilons.  "reference" '\
             'scales all objective values relative to the '\
             'reference point.  This mode is invalid if there '\
             'are any zeros in the reference point.  The '\
             '"reference" mode is compatible with the --integer '\
             'flag, as long as the reference point does not fall '\
             'in a zero-index epsilon box for any objective.')

    # delimiters
    delimiters = parser.add_mutually_exclusive_group()
    delimiters.add_argument('-d', '--delimiter', type=str, default=' ',
        help='input column delimiter, default to space (" ")')
    delimiters.add_argument('--tabs', action='store_true',
        help='use tabs as delimiters')

    # separating and identifying sets
    parser.add_argument('-s', '--separator', type=str,
        help='String which, if found at the beginning of a line, '\
             'causes that line to be treated as a separator between '\
             'solution sets.  To ignore separators instead, treat '\
             'them as comments using one of the comment options.')
    parser.add_argument('--separator-regex', type=regex,
        help='regular expression, which, if a line matches it, '\
             'causes that line to be treated as a separator between '\
             'solution sets')
    parser.add_argument('--no-separate-files', action='store_true',
        help='Treat all inputs as a continuous stream, rather than as '\
             'separate files.  By default, soluion sets will be '\
             'distinguished by file of origin as well as separators '\
             'and index columns.')
    parser.add_argument('--index-columns', type=intrange, nargs='+',
        help='Indices of columns that identify what solution set a point '\
             'belongs to.  Respects --reverse-column-indices.'\
             'If used in combination with separators, the number of '\
             'separators preceding the point in the input will also '\
             'determine what solution set the point belongs to.')
    parser.add_argument('--index-column-names', type=str, nargs='+',
        help='Names of index columns used to identify what solution set a '\
             'point belongs to.  If used in combination with --index-columns, '\
             'maps those column indices to the specified nams.  If used in '\
             'combination with headers, index '\
             'columns are determined '\
             'instead by '\
             'matching names in the header.  If neither index columns '\
             'nor header row are specified, this option is invalid.')
    parser.add_argument('--no-order', action='store_true',
        help='Do not assume that input rows are grouped by set. '\
             'If index columns '\
             'are specified, they alone will be used to distinguish '\
             'sets within an input file.  Unless --no-separate-files '\
             'is also specified, --no-order will still distinguish sets '\
             'by input file.  '\
             'Use of this argument may result in dramatically '\
             'larger memory requirements as well as performance degradation.')

    # comments
    parser.add_argument("-c", "--comment", type=str,
       help="ignore lines starting with this character")
    parser.add_argument("--comment-regex", type=regex,
        help='regular expression, which, if a line matches it, '\
             'causes that line to be treated as a comment and ignored.')

    # Stuff at the top of the file
    parser.add_argument('--skip-initial-lines', type=int, default=0,
        help='Number of lines to skip at the beginning of the file.')
    parser.add_argument('--header-line', type=int,
        help='Input line (counting from 1) on which to find header '\
             'information for the points in the solution sets. '\
             'Implies --skip-initial-lines up to the one on which '\
             'the header is found.  Compatible with a range for '\
             '--skip-initial-lines that includes the header.')
    parser.add_argument('--header-leading-characters', type=str,
        help='Leading characters that identify a new header and '\
             'which are removed, including any whitespace that '\
             'immediately follows, before the header is split '\
             'by a delimiter.')

    # Handling defective input
    parser.add_argument('--malformed-lines', default='warn',
        choices=['ignore', 'empty', 'warn', 'exception'],
        help='How to handle sets with a point that has objectives that '\
             'either cannot be found or cannot be interpreted as numbers. '\
             '"ignore" means skipping malformed input lines and proceeding '\
             'as if nothing were wrong. "empty" means reporting hypervolume '\
             'for a set with malformed input lines as if the set had no '\
             'solutions in it.  "exception" means raising an exception and '\
             'terminating hypervolume computation immediately.' \
             '"warn" behaves like "empty", but also prints a '\
             'message to stderr indicating the defective '\
             'input line and set.')
    parser.add_argument('--empty-set-hypervolume', default='skip',
        choices=['nan', 'zero', 'skip', 'skip-noincrement'],
        help='How to report the hypervolume of a solution set that '\
             'has no solutions in it or is treated as such due to '\
             'malformed input. '\
             '"nan" means that the output will report "NaN" '\
             "for the empty solution set's hypervolume.  "\
             '"zero" means that the output will report 0 '\
             "for the empty solution set's hypervolume.  "\
             '"skip" means that the output will increment the set counter '\
             'if appropriate, but will not contain a row reporting '\
             'hypervolume for the empty solution set.  '\
             '"skip-noincrement" means that the empty solution set '\
             'will be treated as if it did not exist at all.  The set '\
             'counter will not be incremented even if it would '\
             'be appropriate to do so.')

    parser.add_argument("--ZN", action='store_true',
            help='use ZN hypervolume algorithm rather than WFG')
    args = parser.parse_args(argv)
    args.objectives = rerange(args.objectives)
    args.maximize = rerange(args.maximize)
    args.index_columns = rerange(args.index_columns)

    if args.reverse_column_indices:
        if args.objectives is not None:
            args.objectives = [-1 - ob for ob in args.objectives]
        if args.maximize is not None:
            args.maximize = [-1 -ob for ob in args.maximize]
        if args.index_columns is not None:
            args.index_columns = [-1 - c for c in args.index_columns]

    # maximize columns must be objectives
    if args.maximize is not None and args.objectives is not None:
        for index in args.maximize:
            if index not in args.objectives:
                args.objectives.append(index)
    if args.maximize_column_names is not None and \
       args.objective_column_names is not None:
        for name in args.maximize_column_names:
            if name not in args.objective_column_names:
                args.objective_column_names.append(name)

    if args.tabs:
        args.delimiter = "\t"

    return args

def _pure_lines_from_files(files, **kwargs):
    """
    Implementation of lines_from_files that assumes every line is 
    pure data
    """
    header = None
    sep = 0
    name = None
    number = 0
    for fp in files:
        name = fp.name
        number = 0
        header = None
        about = {'sep': sep, 'name': name, 'header': header}
        for line in fp:
            number += 1
            yield (line, number, about)

def _noregex_lines_from_files(files, **kwargs):
    """
    Implementation of lines_from_files that leaves out regex matching.
    kwargs:
    separator
    skip_initial_lines
    comment
    header_line
    header_leading_characters
    """
    header = None
    sep = 0
    name = None
    number = 0
    separator = kwargs.get('separator', None)
    comment = kwargs.get('comment', None)
    if kwargs.get('header_leading_characters', None) is not None:
        header_leader = kwargs['header_leading_characters']
    else:
        header_leader = None

    for fp in files:
        header = None
        sep = 0
        name = fp.name
        number = 0

        if kwargs.get('header_line', None) is not None:
            # lines start with 1 but columns with 0.  sorry!
            for _ in range(1, kwargs.get('header_line')):
                number += 1
                next(fp)
            number += 1
            header = next(fp)
        if kwargs.get('skip_initial_lines', None) is not None:
            for _ in range(number, kwargs.get('skip_initial_lines')):
                number += 1
                next(fp)
        if header_leader is not None and header is not None:
            if header.startswith(header_leader):
                header = header.partition(header_leader)[2].strip()

        about = {'sep': sep, 'name': name, 'header': header}
        for line in fp:
            number += 1
            if header_leader is not None and line.startswith(header_leader):
                sep += 1
                header = line.partition(header_leader)[2].strip()
                about = {'sep': sep, 'name': name, 'header': header}
            elif separator is not None and line.startswith(separator):
                sep += 1
                about = {'sep': sep, 'name': name, 'header': header}
            elif comment is not None and line.startswith(comment):
                pass 
            else:
                yield (line, number, about)

def _regex_lines_from_files(files, **kwargs):
    """
    Implementation of lines_from_files that includes regex matching.
    Definitely slower.
    kwargs:
    separator
    separator_regex
    skip_initial_lines
    comment
    comment_regex
    header_line
    """
    # just pull lines from the noregex version and alter them
    # this saves us having to implement the initial line skipping
    # twice
    noregex = _noregex_lines_from_files(files, **kwargs)
    separator = kwargs.get('separator_regex', None)
    comment = kwargs.get('comment_regex', None)
    sep_offset = 0
    sep = 0
    name = None

    for line, number, about in noregex:
        if about['name'] != name: # next file
            name = about['name']
            sep_offset = 0
            sep = about['sep']
            new_header = about['header']
            new_about = {'sep': sep + sep_offset,
                         'name': name,
                         'header': new_header}

        if about['sep'] != sep:
            sep = about['sep']
            new_about = {'sep': sep + sep_offset,
                         'name': name,
                         'header': new_header}

        if comment is not None:
            if comment.search(line) is not None:
                continue # skip comment line
        if separator is not None:
            if separator.search(line) is not None:
                sep_offset += 1
                new_about = {'sep': sep + sep_offset,
                             'name': name,
                             'header': new_header}
        yield (line, number, new_about)

def lines_from_files(files, **kwargs):
    """
    Generator function yielding lines from files, 
    decorated with file info, separator count in file, and unparsed header.
    kwargs: (see get_args for descriptions)
    separator
    separator_regex
    skip_initial_lines
    comment
    comment_regex
    header_line
    """
    regex = any(kwargs.get(key, None) is not None for key 
                in ['separator_regex', 'comment_regex'])
    impure = any(kwargs.get(key, None) is not None for key
                 in ['separator', 'skip_initial_lines', 'comment',
                     'header_line'])
    if regex:
        implementation = _regex_lines_from_files
    elif impure:
        implementation = _noregex_lines_from_files
    else:
        implementation = _pure_lines_from_files

    decolines = implementation(files, **kwargs)
    for decoline in decolines:
        yield decoline

def rows_from_lines(decolines, **kwargs):
    """
    produces decorated rows, no conversion to numbers yet
    decolines: decorated lines
    kwargs:
    delimiter
    tabs
    """
    if kwargs.get('tabs', False) is True:
        delimiter = "\t"
    else:
        delimiter = kwargs.get('delimiter', ' ')
    header = None
    headerrow = None
    old_about = None
    for line, number, about in decolines:
        line = line.strip()
        if about['header'] is None and header is None:
            pass # no change to header
        elif about['header'] != header:
            header = about['header']
            headerrow = header.split(delimiter)
        if about != old_about:
            old_about = about
            new_about = {'header': headerrow, 
                         'name': about['name'], 
                         'sep': about['sep']}
        row = line.split(delimiter)
        yield(row, number, new_about)

def _ordered_rowsets_from_rows(decorows, **kwargs):
    """
    split rowsets on file, separator, and index columns
    kwargs:
    no_separate_files: grouping omits file and keeps incrementing sep instead
    index_columns
    index_column_names
    """
    if kwargs.get('no_separate_files', False) is True:
        separate = False
    else:
        separate = True

    index_columns = None
    index_column_names = None
    unknown_index = []
    if kwargs.get('index_columns', None) is not None:
        index_columns = kwargs['index_columns']
        unknown_index = [''] * len(index_columns)
    elif kwargs.get('index_column_names', None) is not None:
        index_column_names = kwargs['index_column_names']
        unknown_index = [''] * len(index_column_names)

    name = ''
    sep = 0
    runningsep = 0
    index = [x for x in unknown_index] #copy
    if index_columns is not None:
        unknown_index_columns = [x for x in index_columns]
    else:
        unknown_index_columns = [None for _ in unknown_index]
    current_index_columns = [x for x in unknown_index_columns]
    header = None
    grouping = {'sep': sep, 'name': name, 'index': index}
    inhibit = True # don't yield the dummy at the start

    rowset = []
    for row, number, about in decorows:
        wrapup = False
        if about['name'] != name:
            if separate is False:
                runningsep += sep
            name = about['name']
            sep = about['sep']
            wrapup = True
        if about['sep'] != sep:
            sep = about['sep']
            wrapup = True
        if index_column_names is not None:
            if about['header'] != header:
                header = about['header']
                current_index_columns = [x for x in unknown_index_columns]
                for i, column_name in enumerate(index_column_names):
                    if header is None:
                        break
                    try:
                        current_index_columns[i] = header.index(column_name)
                    except ValueError:
                        current_index_columns[i] = None
            about['index_columns'] = current_index_columns
        if index_columns is not None:
            about['index_columns'] = current_index_columns
        if current_index_columns is not None:
            current_index = [x for x in unknown_index] # copy
            for i,j in enumerate(current_index_columns):
                try:
                    current_index[i] = row[j]
                except TypeError: #None
                    pass # let it remain the empty string
            if not all(x == y for x, y in zip(index, current_index)):
                index = current_index
                wrapup = True
        if wrapup:
            if separate is False:
                del grouping['name']
            grouping['index'] = tuple(grouping['index'])
            if not inhibit:
                yield (rowset, grouping)
            inhibit = False
            rowset = []
            name = about['name']
            sep = about['sep']
            grouping = {
                'sep': runningsep + sep, 
                'name': name, 
                'index': index}
        rowset.append((row, number, about))
    if separate is False:
        del grouping['name']
    grouping['index'] = tuple(grouping['index'])
    if not inhibit:
        yield (rowset, grouping)

def _no_order_rowsets_from_rows(decorows, **kwargs):
    """
    collect all rowsets, but then group only by index and,
    unless no-separate-files is set, by file
    """
    ordered = _ordered_rowsets_from_rows(decorows, **kwargs)
    together = dict()
    for (rowset, grouping) in ordered:
        if kwargs.get('no_separate_files', True) is False:
            key = (grouping['index'], grouping['name'])
        else:
            key = grouping['index']
        collector = together.get(key, [])
        collector.extend(rowset)
        together[key] = collector
    if kwargs.get('no_separate_files', True) is False:
        return [(collector, {'index': index, 'name': name})
                for (index, name), collector in together.iteritems()]
    else:
        return [(collector, {'index': index})
                for index, collector in together.iteritems()]

def rowsets_from_rows(decorows, **kwargs):
    """
    decorows: string row, number, about, where about is
              a dict with headerrow
    kwargs:
    index_columns
    index_column_names
    no_separate_files
    no_order
    """
    if kwargs.get('no_order', False) is True:
        return _no_order_rowsets_from_rows(decorows, **kwargs)
    return _ordered_rowsets_from_rows(decorows, **kwargs)

def _error_ignore(msg): """ ignore an error """; pass
def _error_warn(msg): 
    """ emit a warning on error """
    sys.stderr.write(msg)
    sys.stderr.write("\n")
def _error_input_error(msg):
    """ raise an exception on error """
    raise InputError(msg)

def _convert_objectives_from_indexed_columns(decorowset, **kwargs):
    """
    decorowset: (rowset, grouping)
    kwargs:
    error_handler: Callable taking one argument, a message.
                   The error handler is called for malformed lines.
    """
    column_indices = kwargs['objectives']
    rowset, grouping = decorowset
    error_handler = kwargs.get('error_handler', _error_input_error)

    if len(rowset) == 0:
        return ([], grouping)
    converted = []
    for row, number, about in rowset:
        try:
            for i in column_indices:
                row[i] = float(row[i])
        except ValueError:
            msg = "Could not convert objective value to float "\
                  "in line {0} of {1}."
            msg = msg.format(number, about['name'])
            error_handler(msg)
            continue
        except IndexError:
            msg = "Could not find all sepcified objectives "\
                  "in line {0} of {1}."
            msg = msg.format(number, about['name'])
            error_handler(msg)
            continue
        converted.append((row, number, about))
    return (converted, grouping)

def _convert_objectives_from_named_columns(decorowset, **kwargs):
    """
    decorowset: (rowset, grouping)
    """
    column_names = kwargs['objective_column_names']
    rowset, grouping = decorowset
    if len(rowset) == 0:
        return ([], grouping)
    # we assume that all of the rows in the set have the
    # same header because we split on header changes
    _, number, about = rowset[0]
    header = about['header']
    try:
        column_indices = [header.index(column_name) 
                          for column_name in column_names]
    except ValueError:
        msg = "Header for line {0} in {1} (and {2} more lines) "\
              "does not contain all of the specified objectives"
        msg = msg.format(number, about['name'], len(rowset) - 1)
        raise InputError(msg)
    return _objectives_from_indexed_columns(
        decorowset, objectives=column_indices, **kwargs)

def _convert_objectives_from_all_columns(decorowset, **kwargs):
    """
    All rows need to have the same length.
    decorowset: (rowset, grouping)
    """
    rowset, grouping = decorowset
    if len(rowset) == 0:
        return ([], grouping)
    error_handler = kwargs.get('error_handler')
    
    converted = []
    firstrow = True
    nobj=0

    if kwargs.get('index_columns', None) is not None:
        index_columns = kwargs['index_columns']
        for row, number, about in rowset:
            nobj_seen = 0
            try:
                for i, x in enumerate(row):
                    if i in index_columns:
                        continue
                    nobj_seen += 1
                    row[i] = float(x)
            except ValueError:
                msg = "Could not convert objective value to float "\
                      "in line {0} of {1}."
                msg = msg.format(number, about['name'])
                error_handler(msg)
                continue
            if firstrow:
                nobj = nobj_seen
                firstrow = False
            if nobj_seen != nobj:
                msg = "Wrong number of objectives "\
                      "in line {0} of {1}."
                msg = msg.format(number, about['name'])
                error_handler(msg)
                continue
            converted.append((row, number, about))
    elif kwargs.get('index_column_names', None) is not None:
        for row, number, about in rowset:
            index_columns = about['index_columns']
            nobj_seen = 0
            try:
                for i, x in enumerate(row):
                    if i in index_columns:
                        continue
                    nobj_seen += 1
                    row[i] = float(x)
            except ValueError:
                msg = "Could not convert objective value to float "\
                      "in line {0} of {1}."
                msg = msg.format(number, about['name'])
                error_handler(msg)
                continue
            if firstrow:
                nobj = nobj_seen
                firstrow = False
            if nobj_seen != nobj:
                msg = "Wrong number of objectives "\
                      "in line {0} of {1}."
                msg = msg.format(number, about['name'])
                error_handler(msg)
                continue
            converted.append((row, number, about))
    else: # no index columns, yay!
        for row, number, about in rowset:
            nobj_seen = 0
            try:
                for i, x in enumerate(row):
                    row[i] = float(x)
                    nobj_seen += 1
            except ValueError:
                msg = "Could not convert objective value to float "\
                      "in line {0} of {1}."
                msg = msg.format(number, about['name'])
                error_handler(msg)
                continue
            if firstrow:
                nobj = nobj_seen
                firstrow = False
            if nobj_seen != nobj:
                msg = "Wrong number of objectives "\
                      "in line {0} of {1}."
                msg = msg.format(number, about['name'])
                error_handler(msg)
                continue
            converted.append((row, number, about))
    return (converted, grouping)

def convert_objectives_from_rowsets(decorowsets, **kwargs):
    """
    convert objectives in the manner requested by kwargs
    kwargs:
    Commandline args describe options in detail
    malformed_lines: What to do if input lines are malformed.
    objective_column_names: Names of objective columns
    objectives: Indices of objective columns
    """

    on_malformed = kwargs.get('malformed_lines', 'empty')
    if on_malformed in ['empty', 'exception']:
        error_handler = _error_input_error
    elif on_malformed == 'warn':
        error_handler = _error_warn
    elif on_malformed == 'ignore':
        error_handler = _error_ignore
    if kwargs.get('objective_column_names', None) is not None:
        implementation = _convert_objectives_from_named_columns
    elif kwargs.get('objectives', None) is not None:
        implementation = _convert_objectives_from_indexed_columns
    else:
        implementation = _convert_objectives_from_all_columns
    for decorowset in decorowsets:
        if on_malformed == 'empty':
            try:
                yield implementation(decorowset, error_handler=error_handler, **kwargs)
            except InputError as ie:
                yield ([], decorowset[1])
        else: # allow exception to propagate
            yield implementation(decorowset, error_handler=error_handler, **kwargs)

def hypervolume(rows, **kwargs):
    """
    kwargs:
    objectives: indices to use. if not specified, use all floats in each row
    objective_column_names: names of objective columns.  requires header
    maximize: indices to maximize
    maximize_all: maximize all objectives
    maximize_column_names: maximize named columns.  requires header
    header: column names
    epsilons: epsilons in order of provided objective indices
    integer: Compute hypervolume based on box corners.  Requires epsilons.
    reference: Reference point.  Default to zero if not specified.
    scale: 'none' or None means not to scale.  'epsilon' divides objective
           values by epsilon (requires epsilons).  'reference' divides 
           objective values by reference point but fails if 
           reference point has zeros.
    """
    rows, _ = apply_maximization(rows, **kwargs)
    reference = kwargs.get('reference', None)
    if reference is not None:
        reference, _ = apply_maximization(
            [reference], **kwargs)
        reference = reference[0]
        reference = [x for x in reference] # copy just in case

    # finally, we can drop the stuff that's not objectives!!!!!
    objectives = kwargs.get('objectives', None)
    objective_column_names = kwargs.get('objective_column_names', None)
    header = kwargs.get('header', None)
    if objectives is not None:
        pass # easy way to exclude the condition where both are not None
    elif objective_column_names is not None and header is not None:
        try:
            objectives = [header.index(name) for name 
                          in objective_column_names]
        except ValueError:
            raise InputError('not all objective names found in header')

    objective_rows = []
    if objectives is None:
        for row in rows:
            floats = [x for x in row if isinstance(x, float)]
            objective_rows.append(floats)
    else:
        try:
            for row in rows:
                objective_rows.append([row[i] for i in objectives])
        except IndexError:
            raise InputError('objective missing from input row')

    sys.stdout.flush()
    # now we have objective rows, and maximization has been fixed up
    # validate nobj
    nobj=None
    for row in objective_rows:
        if nobj is None:
            nobj = len(row)
        if nobj != len(row):
            raise InputError('input row has the wrong number of objectives')

    # sort if epsilons specified
    epsilons = None
    if kwargs.get('epsilons', None) is not None:
        epsilons = kwargs['epsilons']
        objective_rows = pareto.eps_sort(objective_rows, epsilons=epsilons)

    # Make integral, including reference point
    # They and math on them are exact per the IEEE standard as long as
    # we stay in integer-land and don't exceed certain limits.
    # See http://stackoverflow.com/questions/3387655/safest-way-to-convert-float-to-integer-in-python
    if kwargs.get('integer', False) is True and epsilons is not None:
        if reference is not None:
            for i, val in enumerate(reference):
                reference[i] = math.floor(val / epsilons[i])
        for row in objective_rows:
            for i, val in enumerate(row):
                row[i] = math.floor(val / epsilons[i])
    
    # apply scale
    scale = kwargs.get('scale', None)
    if scale is not None and scale != 'none':
        if scale == 'epsilon' and epsilons is not None:
            if integer is False: # if true, scale is already applied
                # scale reference point
                if reference is not None:
                    for i, val in enumerate(reference):
                        reference[i] = val / epsilons[i]
                for row in objective_rows:
                    for i, val in enumerate(row):
                        row[i] = val / epsilons[i]
        if scale == 'reference' and reference is not None:            
            if 0.0 in reference:
                raise InputError('cannot scale to zero reference')
            for row in objective_rows:
                for i, val in enumerate(row):
                    row[i] = val / reference[i]
            reference = [1.0] * len(reference)

    # transform to reference point, clamping
    # wfg assumes maximization, so this assures nonnegative values
    if reference is not None:
        if len(reference) == 1:
            reference = reference * len(objective_rows[0])
        elif len(reference) != len(objective_rows[0]):
            raise ReferencePointError('Reference point does not match dimension of set')
        for row in objective_rows:
            for i, val in enumerate(row):
                if val > reference[i]: #clamp!
                    val = reference[i]
                row[i] = reference[i] - val

    # finally!  compute hypervolume
    if kwargs.get('ZN', False) is True:
        hv = zn.zn(objective_rows)
    else:
        hv = wfg.wfg(objective_rows)
    return hv

def nadir(rows):
    """
    Assume each row has the same number of floats in it.
    Assume minimization.
    rows (list of lists): data to convert
    returns nadir point
    """
    nadir = [x for x in rows[0] if isinstance(x, float)]
    for row in rows[1:]:
        current = [x for x in row if isinstance(x, float)]
        nadir = [max([n, x]) for n,x in zip(nadir, current)]
    return nadir

def zenith(rows):
    """
    Assume each row has the same number of floats in it.
    Assume minimization.
    rows (list of lists): data to convert
    returns zenith point
    """
    zenith = [x for x in rows[0] if isinstance(x, float)]
    for row in rows[1:]:
        current = [x for x in row if isinstance(x, float)]
        zenith = [min([n, x]) for n,x in zip(zenith, current)]
    return zenith

def apply_maximization(rows, **kwargs):
    """
    apply maximization and produce kwargs that don't mention it
    """
    maximize = kwargs.get('maximize', None)
    maximize_all = kwargs.get('maximize_all', None)
    maximize_column_names =  kwargs.get('maximize_column_names', None)

    if maximize is None and maximize_column_names is not None:
        if header is not None:
            try:
                maximize = [header.index(col) for col in maximize_column_names]
            except ValueError:
                msg = "Not all maximization columns could be found in header"
                raise InputError(msg)

    flipped_columns = []
    transformed = []
    if maximize is not None:
        for row in rows:
            trow = []
            for i, x in enumerate(row):
                if i in maximize:
                    x = -1.0 * x
                    # which objective in the sequence of objectives was maximized
                    flipped_columns.append(maximize.index(i))
                trow.append(x)
            transformed.append(trow)
    elif maximize_all is True:
        for row in rows:
            trow = []
            for i, x in enumerate(row):
                flipped_columns.append(i)
                if isinstance(x, float):
                    trow.append(-1.0 * x)
                else:
                    trow.append(x)
            transformed.append(trow)
    else: # nothing to maximize, return the input
        transformed = rows
    transformed_kwargs = dict()
    for key, val in kwargs.iteritems():
        if "maximize" not in key:
            transformed_kwargs[key] = val
    transformed_kwargs['flipped_columns'] = flipped_columns
    return transformed, transformed_kwargs

def hypervolumes_from_converted_sets(decosets, **kwargs):
    """
    Float conversion has been done at this point.
    decosets: (rowset, grouping)
    """
    reference =  kwargs.get('reference', None)
    auto_reference = kwargs.get('auto_reference', None)
    if reference is None and auto_reference is None:
        auto_reference = 'zero'

    refpoint = None
    antirefpoint = None
    for rowset, grouping in decosets:
        if len(rowset) == 0:
            yield ('empty', grouping)
            continue
        header = rowset[0][2]['header']
        rows = [row for row, _, _ in rowset]
        # account for maximization now, even though hypervolume also
        # does it
        rows, processed_keywords = apply_maximization(
            rows, header=header, **kwargs)
        
        if reference is not None:
            flipped_columns = processed_keywords['flipped_columns']
            maxreference = []
            for i, x in enumerate(reference):
                if i in flipped_columns:
                    maxreference.append(-1.0 * reference[i])
                else:
                    maxreference.append(reference[i])
        else:
            maxreference = reference

        if auto_reference is None:
            refpoint = maxreference
        elif auto_reference == 'zero':
            refpoint = [0 for x in rows[0] if isinstance(x, float)]
            if reference is not None:
                refpoint = nadir([refpoint, maxreference])
        elif auto_reference == 'first':
            auto_reference = 'fixed'
            refpoint = nadir(rows)
            if reference is not None:
                refpoint = nadir([refpoint, maxreference])
        elif auto_reference == 'fixed':
            pass # just don't update refpoint
        elif auto_reference == 'full':
            if refpoint is not None:
                refpoint = nadir([nadir(rows), refpoint])
            else:
                refpoint = nadir(rows)
            if reference is not None:
                refpoint = nadir([refpoint, maxreference])
            if antirefpoint is not None:
                antirefpoint = zenith([zenith(rows), antirefpoint])
            else:
                antirefpoint = zenith(rows)
            grouping['zenith'] = antirefpoint
        else:
            refpoint = maxreference
        processed_keywords['reference'] = refpoint
        grouping['reference'] = refpoint
        processed_keywords['header'] = header
        try:
            hv = hypervolume(rows, **processed_keywords)
            yield (hv, grouping)
        except ReferencePointError as rpe:
            if reference is not None:
                sys.stderr.write("Supplied reference point does not match all input sets.\n")
            if auto_reference == 'zero':
                sys.stderr.write("THIS SHOULD NOT HAPPEN:\n"\
                                 "An automatic zero reference point should never be the wrong size.\n")
            else:
                sys.stderr.write("WARNING: Sets changed size while using automatic reference point.\n"\
                                 "         Skipping set index '{0}' from filename '{1}' "\
                                 "after separator '{2}'.\n".format(
                                     grouping.get('index', ''), grouping.get('name', ""),
                                     grouping.get('sep', "")))

def write_output(hypervolumes, **kwargs):
    """
    write output
    kwargs:
    no_separate_files
    no_order
    index_column_names
    index_columns
    reference
    auto_reference
    """
    grouping_bits = []
    header = []
    if kwargs.get('no_separate_files', False) is False:
        grouping_bits.append('name')
        header.append('input')
    if kwargs.get('no_order', False) is False:
        grouping_bits.append('sep')
        header.append('front number')
    if kwargs.get('index_column_names', None) is not None:
        grouping_bits.append('index')
        header.extend(kwargs['index_column_names'])
    elif kwargs.get('index_columns', None) is not None:
        grouping_bits.append('index')
        header.extend(["i{0}".format(x) for x in range(len(kwargs.get('index_columns')))])
    if kwargs.get('reference') is not None:
        grouping_bits.append('reference')
        header.append('reference')
    elif kwargs.get('auto_reference') == 'full':
        grouping_bits.append('reference')
        grouping_bits.append('zenith')
        header.append('reference')
        header.append('zenith')
    elif kwargs.get('auto_reference') == 'first':
        grouping_bits.append('reference')
        header.append('reference')
    header.append('hv')
    delimiter = kwargs.get('delimiter')
    output = kwargs.get('output')
    output.write(delimiter.join(header))
    output.write('\n')

    #    choices=['nan', 'zero', 'skip', 'skip-noincrement'],
    on_empty = kwargs.get("empty_set_hypervolume", "nan")
    sep_offset = 0
    current_file = None
    for hv, grouping in hypervolumes:
        effective_grouping = dict([(k, v) for k, v in grouping.items()])
        if on_empty == 'skip-noincrement': 
            name = effective_grouping.get('name', None)
            if name != current_file and kwargs.get('no_separate_files', False) is False:
                current_file = name
                sep_offset = 0
            if effective_grouping.get('sep', None) is not None:
                effective_grouping['sep'] -= sep_offset
        if hv == 'empty':
            if on_empty == 'nan':
                the_hypervolume = float("NaN")
            elif on_empty == 'zero': 
                the_hypervolume = 0.0
            elif on_empty == 'skip':
                continue
            elif on_empty == 'skip-noincrement':
                sep_offset += 1
                continue
        else:
            the_hypervolume=hv # why does this prevent a weird bug?
        output.write(delimiter.join(
            [str(effective_grouping[bit]) for bit in grouping_bits]))
        output.write(delimiter)
        output.write("{0:.5g}".format(the_hypervolume))
        output.write('\n')

def cli(args):
    """
    args (namespace): the result of parsing input options
    """
    files = args.inputs
    kwargs = args.__dict__

    # pipeline: everything gets decorated with identifying info
    # pipeline stage 0: emit lines from files, decorated with relevant info
    lff = lines_from_files(files, **kwargs)
    # pipeline stage 1: process delimiters to produce rows
    rfl = rows_from_lines(lff, **kwargs)
    # pipeline stage 2: collect rows into row sets
    rfr = rowsets_from_rows(rfl, **kwargs)
    # pipeline stage 3: extract objectives from rows
    cfr = convert_objectives_from_rowsets(rfr, **kwargs)

    # pipeline stage 4: compute hypervolume
    hfc = hypervolumes_from_converted_sets(cfr, **kwargs)
    # pipeline stage 5: write outputs
    write_output(hfc, **kwargs)

if __name__ == '__main__':
    cli(get_args(sys.argv))

