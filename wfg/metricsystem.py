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
    parser.add_argument('-m', '--maximize', type=intrange, nargs='+',
        help='Objective columns to maximize, default is none. '\
             'Maximized objectives are multiplied by -1.0 and '\
             'treated as minimization objectives.'\
             'The final hypervolume is also corrected by negation '\
             'if there is an '\
             'odd number of maximization objectives.')
    parser.add_argument('-M', '--maximize-all', action='store_true',
        help='maximize all objectives')
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
             'number of objectives will be evaluated.  Default '\
             'behavior is the same as "-R 0.0", treating the '\
             'origin as a reference point.  Hard-to-interpret '\
             'hypervolume data will result if the data straddle '\
             'the reference point, so a zenith or nadir point is '\
             'recommended.') 

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
    parser.add_argument('-s', '--separator', type=str, default='#',
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
    parser.add_argument('--index-columns', type=intrange,
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
    parser.add_argument('--header-line', type=int, default=-1,
        help='Input line (counting from 0) on which to find header '\
             'information for the points in the solution sets. '\
             'Implies --skip-initial-lines up to the one on which '\
             'the header is found.  Compatible with a range for '\
             '--skip-initial-lines that includes the header.')
    parser.add_argument('--header-regex', type=regex,
        help='Lines matching this regular expression will be '\
             'treated as headers for all following lines, '\
             'until another header line or the end of the file '\
             'is reached.  A new header is considered to separate '\
             'solution sets.')

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
             'if appropriate, but will not contain a row reporting hypervolume '\
             'for the empty solution set.  '\
             '"skip-noincrement" means that the empty solution set '\
             'will be treated as if it did not exist at all.  The set '\
             'counter will not be incremented even if it would '\
             'be appropriate to do so.')

    args = parser.parse_args(argv)
    args.objectives = rerange(args.objectives)
    args.maximize = rerange(args.maximize)

    if args.reverse_column_indices:
        if args.objectives is not None:
            args.objectives = [-1 - ob for ob in args.objectives]
        if args.maximize is not None:
            args.maximize = [-1 -ob for ob in args.maximize]

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
    """
    header = None
    sep = 0
    name = None
    number = 0
    separator = kwargs.get('separator', None)
    comment = kwargs.get('comment', None)

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

        about = {'sep': sep, 'name': name, 'header': header}
        if separator is None and comment is None:
            for line in fp:
                number += 1
                yield (line, number, about)
        elif separator is None:
            for line in fp:
                number += 1
                if not line.startswith(comment):
                    yield (line, number, about)
        elif comment is None:
            for line in fp:
                number += 1
                if line.startswith(separator):
                    sep += 1
                    about = {'sep': sep, 'name': name, 'header': header}
                else:
                    yield (line, number, about)
        else:
            for line in fp:
                number += 1
                if not line.startswith(comment):
                    if line.startswith(separator):
                        sep += 1
                        about = {'sep': sep, 'name': name, 'header': header}
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
    header_regex
    """
    # just pull lines from the noregex version and alter them
    # this saves us having to implement the initial line skipping
    # twice
    noregex = _noregex_lines_from_files(files, **kwargs)
    separator = kwargs.get('separator_regex', None)
    comment = kwargs.get('comment_regex', None)
    header = kwargs.get('header_regex', None)
    sep_offset = 0
    name = None
    override_header = None

    for line, number, about in noregex:
        if about.name != name: # next file
            name = about.name
            sep_offset = 0
            override_header = None
        if comment is not None:
            if comment.search(line) is not None:
                continue # skip comment line
        if separator is not None:
            if separator.search(line) is not None:
                sep_offset += 1
        if header is not None:
            if header.search(line) is not None:
                override_header = line
                sep_offset += 1
                continue
        about['sep'] = about['sep'] + sep_offset
        if override_header is not None:
            about['header'] = override_header
        yield (line, number, about)

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
    header_regex
    """
    regex = any(kwargs.get(key, None) is not None for key 
                in ['separator_regex', 'comment_regex', 'header_regex'])
    impure = any(kwargs.get(key, None) is not None for key
                 in ['separator', 'skip_initial_lines', 'comment',
                     'header_line'])
    print('separator {0}'.format(kwargs.get('separator')))
    print('comment {0}'.format(kwargs.get('comment')))
    print('header_line {0}'.format(kwargs.get('header_line')))
    if regex:
        print('using regex')
        implementation = _regex_lines_from_files
    elif impure:
        print('using impure')
        implementation = _noregex_lines_from_files
    else:
        print('using pure')
        implementation = _pure_lines_from_files

    lines = implementation(files, **kwargs)
    for line in lines:
        yield line

def cli(args):
    """
    args (namespace): the result of parsing input options
    """
    # pipeline stage 0: emit lines from files, decorated with relevant info
    files = args.inputs
    kwargs = args.__dict__
    lff = lines_from_files(files, **kwargs)
    for line, number, about in lff:
        print(line, number, about) 
    # pipeline stage 1: emit line sets decorated with identifying info
    # pipeline stage 2: emit row sets (objectives only, converted to float) 
    #                   decorated with identifying info
    # pipeline stage 2a: emit eps-nondom sorted row sets, same decoration
    # pipeline stage 2b: emit integral row sets, same decoration.
    # pipeline stage 3: emit hypervolumes
    # pipeline stage 4: write outputs

if __name__ == '__main__':
    cli(get_args(sys.argv))

