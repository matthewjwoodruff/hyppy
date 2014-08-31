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

def get_args(argv):
    """ Get command-line arguments """
    prog = argv.pop(0)
    parser = argparse.ArgumentParser(
        prog=prog, description='Metrics Calculation for Pareto Approximations')
    parser.add_argument('inputs', type=argparse.FileType('r'), nargs='+',
        help='input filenames, use - for standard input')
    parser.add_argument('-o', '--objectives', type=intrange, nargs='+',
        help='objective columns (zero-indexed)')
    parser.add_argument('-m', '--maximize', type=intrange, nargs='+',
        help='objective columns to maximize')
    parser.add_argument('-M', '--maximize-all', action='store_true',
        help='maximize all objectives')
    parser.add_argument('--output', type=argparse.FileType('w'),
        default=sys.stdout,
        help='output filename, defaulting to stdout')

    delimiters = parser.add_mutually_exclusive_group()
    delimiters.add_argument('-d', '--delimiter', type=str, default=' ',
        help='input column delimiter, default to space (" ")')
    delimiters.add_argument('--tabs', action='store_true',
        help='use tabs as delimiters')

    # what else we need
    # * something to separate sets
    # * something to identify sets
    # * header, comments

    return parser.parse_args(argv)
