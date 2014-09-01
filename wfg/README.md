# metricsystem.py

```
usage: metricsystem.py [-h] [--output OUTPUT] [-o OBJECTIVES [OBJECTIVES ...]]
                       [--objective-column-names OBJECTIVE_COLUMN_NAMES [OBJECTIVE_COLUMN_NAMES ...]]
                       [-m MAXIMIZE [MAXIMIZE ...]] [-M]
                       [--maximize-column-names MAXIMIZE_COLUMN_NAMES [MAXIMIZE_COLUMN_NAMES ...]]
                       [--reverse-column-indices] [-e EPSILONS [EPSILONS ...]]
                       [--integer] [-R REFERENCE [REFERENCE ...]]
                       [--auto-reference {zero,first,full}]
                       [--scale {none,epsilon,reference}]
                       [-d DELIMITER | --tabs] [-s SEPARATOR]
                       [--separator-regex SEPARATOR_REGEX]
                       [--no-separate-files]
                       [--index-columns INDEX_COLUMNS [INDEX_COLUMNS ...]]
                       [--index-column-names INDEX_COLUMN_NAMES [INDEX_COLUMN_NAMES ...]]
                       [--no-order] [-c COMMENT]
                       [--comment-regex COMMENT_REGEX]
                       [--skip-initial-lines SKIP_INITIAL_LINES]
                       [--header-line HEADER_LINE]
                       [--header-leading-characters HEADER_LEADING_CHARACTERS]
                       [--malformed-lines {ignore,empty,warn,exception}]
                       [--empty-set-hypervolume {nan,zero,skip,skip-noincrement}]
                       inputs [inputs ...]

Metrics Calculation for Pareto Approximations

positional arguments:
  inputs                input filenames, use - for standard input

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT       output filename, defaulting to stdout
  -o OBJECTIVES [OBJECTIVES ...], --objectives OBJECTIVES [OBJECTIVES ...]
                        objective columns (zero-indexed), default is all
  --objective-column-names OBJECTIVE_COLUMN_NAMES [OBJECTIVE_COLUMN_NAMES ...]
                        Names of objective columns
  -m MAXIMIZE [MAXIMIZE ...], --maximize MAXIMIZE [MAXIMIZE ...]
                        Columns to maximize, default is none. All columns for
                        which maximization is specified are considered to be
                        objectives. Maximized objectives are multiplied by
                        -1.0 and treated as minimization objectives.The final
                        hypervolume is also corrected by negation if there is
                        an odd number of maximization objectives.
  -M, --maximize-all    maximize all objectives
  --maximize-column-names MAXIMIZE_COLUMN_NAMES [MAXIMIZE_COLUMN_NAMES ...]
                        Names of columns to maximize. These are automatically
                        considered to be objectives
  --reverse-column-indices
                        Reverse the order of column indices. May be useful if
                        your objectives are at the end of a row of unknown
                        length. Make sure -e and -m are consistent with the
                        order you specify.
  -e EPSILONS [EPSILONS ...], --epsilons EPSILONS [EPSILONS ...]
                        Epsilons, one per objective. Specifying epsilons
                        causes data to be passed through an epsilon-
                        nondomination sort before hypervolume is computed.
  --integer             Use an integer version of the hypervolume computation,
                        based on epsilon box index relative to reference
                        point's epsilon box index. Epsilon box indexes are
                        always referenced to 0 so that box boundaries are
                        consistent with those used in the sort. Specifying ths
                        option will significantly change the computed
                        hypervolumes. This option has no effect if epsilons
                        are not specified.
  -R REFERENCE [REFERENCE ...], --reference REFERENCE [REFERENCE ...]
                        Reference point. If one value is specified, it is
                        duplicated for all objectives. If multiple values are
                        specified, only solution sets with the same number of
                        objectives will be evaluated. Reference point must be
                        a nadir. Any inferior points will be clamped to the
                        reference.
  --auto-reference {zero,first,full}
                        Automatic reference point behavior. "zero" means that
                        the origin is used. "first" means that the nadir point
                        for the first solution set in the input is used as a
                        reference for every solution set. "full" means that
                        the nadir point for all sets combined is used to
                        compute hypervolume for each. The full option means
                        that output will be written twice, so that an aborted
                        run still produces useable output. The first-pass
                        output will include the nadir point used for each
                        set.Combining auto reference with a specified
                        reference point means that the reference point will be
                        the nadir of the specified and automatic points.If
                        neither this option nor a reference point is
                        specified, the "zero" behavior is applied.
  --scale {none,epsilon,reference}
                        Scaling mode. "none" is the default. "epsilon" is
                        synonymous with "none" unless epsilons are specified,
                        in which case all objective values are scaled by their
                        epsilons. "reference" scales all objective values
                        relative to the reference point. This mode is invalid
                        if there are any zeros in the reference point. The
                        "reference" mode is compatible with the --integer
                        flag, as long as the reference point does not fall in
                        a zero-index epsilon box for any objective.
  -d DELIMITER, --delimiter DELIMITER
                        input column delimiter, default to space (" ")
  --tabs                use tabs as delimiters
  -s SEPARATOR, --separator SEPARATOR
                        String which, if found at the beginning of a line,
                        causes that line to be treated as a separator between
                        solution sets. To ignore separators instead, treat
                        them as comments using one of the comment options.
  --separator-regex SEPARATOR_REGEX
                        regular expression, which, if a line matches it,
                        causes that line to be treated as a separator between
                        solution sets
  --no-separate-files   Treat all inputs as a continuous stream, rather than
                        as separate files. By default, soluion sets will be
                        distinguished by file of origin as well as separators
                        and index columns.
  --index-columns INDEX_COLUMNS [INDEX_COLUMNS ...]
                        Indices of columns that identify what solution set a
                        point belongs to. Respects --reverse-column-indices.If
                        used in combination with separators, the number of
                        separators preceding the point in the input will also
                        determine what solution set the point belongs to.
  --index-column-names INDEX_COLUMN_NAMES [INDEX_COLUMN_NAMES ...]
                        Names of index columns used to identify what solution
                        set a point belongs to. If used in combination with
                        --index-columns, maps those column indices to the
                        specified nams. If used in combination with headers,
                        index columns are determined instead by matching names
                        in the header. If neither index columns nor header row
                        are specified, this option is invalid.
  --no-order            Do not assume that input rows are grouped by set. If
                        index columns are specified, they alone will be used
                        to distinguish sets within an input file. Unless --no-
                        separate-files is also specified, --no-order will
                        still distinguish sets by input file. Use of this
                        argument may result in dramatically larger memory
                        requirements as well as performance degradation.
  -c COMMENT, --comment COMMENT
                        ignore lines starting with this character
  --comment-regex COMMENT_REGEX
                        regular expression, which, if a line matches it,
                        causes that line to be treated as a comment and
                        ignored.
  --skip-initial-lines SKIP_INITIAL_LINES
                        Number of lines to skip at the beginning of the file.
  --header-line HEADER_LINE
                        Input line (counting from 1) on which to find header
                        information for the points in the solution sets.
                        Implies --skip-initial-lines up to the one on which
                        the header is found. Compatible with a range for
                        --skip-initial-lines that includes the header.
  --header-leading-characters HEADER_LEADING_CHARACTERS
                        Leading characters that identify a new header and
                        which are removed, including any whitespace that
                        immediately follows, before the header is split by a
                        delimiter.
  --malformed-lines {ignore,empty,warn,exception}
                        How to handle sets with a point that has objectives
                        that either cannot be found or cannot be interpreted
                        as numbers. "ignore" means skipping malformed input
                        lines and proceeding as if nothing were wrong. "empty"
                        means reporting hypervolume for a set with malformed
                        input lines as if the set had no solutions in it.
                        "exception" means raising an exception and terminating
                        hypervolume computation immediately."warn" behaves
                        like "empty", but also prints a message to stderr
                        indicating the defective input line and set.
  --empty-set-hypervolume {nan,zero,skip,skip-noincrement}
                        How to report the hypervolume of a solution set that
                        has no solutions in it or is treated as such due to
                        malformed input. "nan" means that the output will
                        report "NaN" for the empty solution set's hypervolume.
                        "zero" means that the output will report 0 for the
                        empty solution set's hypervolume. "skip" means that
                        the output will increment the set counter if
                        appropriate, but will not contain a row reporting
                        hypervolume for the empty solution set. "skip-
                        noincrement" means that the empty solution set will be
                        treated as if it did not exist at all. The set counter
                        will not be incremented even if it would be
                        appropriate to do so.
```
