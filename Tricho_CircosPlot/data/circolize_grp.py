#!/usr/bin/env python


# Copyright (C) 2017  Shengwei Hou, housw2010'at'gmail'dot'com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function
import argparse
import os
import sys
import numpy as np


def _get_headers(path2grp):
    headers = []
    with open(path2grp, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                headers.append(line.strip())
            else:
                break
    return headers


def _get_ncol(headers):
    ncol = None
    for line in headers:
        line = line.lstrip("# ")
        if not ncol:
            ncol = len(line.strip().split())-1
        else:
            assert ncol == len(line.strip().split())-1, "lines in headers have different number of columns!"
    return ncol


def get_coverage_from_grp(path2grp):
    """
    :param path2grp: input grp file
    :return:         the headers and a list of grp numpy arrays.

    This function used to read in grp file, the common 2 column input grp file format like this:
        # BASE fwd_grp rev_grp
        # colour 0:204:0 0:0:255
        1 0.00 -0.00
        2 0.00 -0.00
        3 0.00 -0.00

    or the 4 column grp file format like this:
        # BASE fwd_cov fwd_tss rev_cov rev_tss
        # colour 0:204:0 255:0:0 0:0:255 255:0:0
        1 0.00 0.00 -0.00 -0.00
        2 0.00 0.00 -0.00 -0.00
        3 0.00 0.00 -0.00 -0.00
    """
    headers = _get_headers(path2grp)
    ncol = _get_ncol(headers)

    fwd_cov = []
    rev_cov = []

    with open(path2grp, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            if ncol == 2:
                fwd_cov.append(abs(float(line[1])))
                rev_cov.append(abs(float(line[2])))
            else:
                assert ncol == 4, "input grp file should have either 2 or 4 data columns!"
                fwd_cov.append(abs(float(line[1])))
                rev_cov.append(abs(float(line[3])))

    return (np.array(fwd_cov), np.array(rev_cov))


def grp2circos(input_grp, chromosome_name, out_plus_coverage, out_minus_coverage, window_size=1000):
    """ given an input grp file, calculate average coverages on the plus and minus strands,
        write out to out_plus_coverage and out_minus_coverage in bed format
    """

    with open(out_plus_coverage, "w") as oh1, open(out_minus_coverage, "w") as oh2:
        oh1.write("#chromosome\tstart\tend\tcoverage\n")
        oh2.write("#chromosome\tstart\tend\tcoverage\n")

        fwd_cov, rev_cov = get_coverage_from_grp(input_grp)
        seq_length = len(fwd_cov)
        name  = chromosome_name

        for window_start in range(0, seq_length, window_size):
            window_end = window_start + window_size \
                         if window_start + window_size < seq_length \
                         else seq_length
            fwd_window_cov = fwd_cov[window_start:window_end].mean()
            rev_window_cov = rev_cov[window_start:window_end].mean()
            oh1.write( name + "\t"+ str(window_start) +"\t"+ str(window_end) +"\t" +str(fwd_window_cov) +"\n")
            oh2.write( name + "\t"+ str(window_start) +"\t"+ str(window_end) +"\t" +str(rev_window_cov) +"\n")


def main():


    # parser
    parser = argparse.ArgumentParser(description="create plus and minus strand\
                         coverage data for circos visualization")
    parser.add_argument("input_file", help="input grp file")
    parser.add_argument("chromosome_name", help="chromosome name for input grp file")
    parser.add_argument("-w", "--window_size", type=int, default=1000,
                        help="window size to calculate averaged coverage" \
                        "default=1000")
    parser.add_argument("-p", "--prefix", help="output prefix")
    parser.add_argument("-o", "--out_folder", help="output directory, default=./", default="./")
    parser.add_argument("-f", "--force", action="store_true", help="force to overwrite the output file")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    if len(sys.argv) < 2:
        print("\nERROR: Not enough parameters were provided, please refer to the usage.\n", file=sys.stderr)
        print(parser.format_help(), file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # input and output handeling
    if args.prefix:
        prefix = args.prefix
    else:
        basename = os.path.basename(args.input_file)
        prefix = os.path.splitext(basename)[0]
    out_plus_coverage = os.path.join(args.out_folder, prefix+"_plus_coverage.txt")
    out_minus_coverage = os.path.join(args.out_folder, prefix+"_minus_coverage.txt")

    if os.path.exists(out_plus_coverage) or os.path.exists(out_minus_coverage):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # convert
    grp2circos(args.input_file, args.chromosome_name,
               out_plus_coverage, out_minus_coverage,
               args.window_size)


if __name__ == "__main__":
    main()
