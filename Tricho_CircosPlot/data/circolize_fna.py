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
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC_skew


def fna2circos(input_fna, out_GC_content, out_GC_skew, window_size=1000):
    """ given an input fna file, calculate GC content and GC skew using
        BioPython at defined window_size, write out to out_GC_content and
        out_GC_skew files
    """
    with open(out_GC_content, "w") as oh1, open(out_GC_skew, "w") as oh2:
        oh1.write("#chromosome\tstart\tend\tgc_content\n")
        oh2.write("#chromosome\tstart\tend\tgc_skew\n")
        for record in SeqIO.parse(input_fna, "fasta"):
            name = record.name
            seq = record.seq
            seq_length = len(seq)
            for window_start in range(0, len(seq), window_size):
                window_end = window_start + window_size \
                             if window_start + window_size < seq_length \
                             else seq_length
                window_seq = seq[window_start:window_end]
                gc_content = GC(window_seq)
                gc_skew = GC_skew(window_seq, window=len(window_seq))[0]
                oh1.write( name + "\t"+ str(window_start) +"\t"+ str(window_end) +"\t" +str(gc_content) +"\n")
                oh2.write( name + "\t"+ str(window_start) +"\t"+ str(window_end) +"\t" +str(gc_skew) +"\n")


def main():


    # parser
    parser = argparse.ArgumentParser(description="create GC_content and GC_skew \
                         data for circos visualization using Biopython")
    parser.add_argument("input_file", help="input fna file")
    parser.add_argument("-w", "--window_size", type=int, default=1000,
                        help="window size to calculate GC content and GC skew" \
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
    out_GC_content = os.path.join(args.out_folder, prefix+"_GC_content.txt")
    out_GC_skew = os.path.join(args.out_folder, prefix+"_GC_skew.txt")

    if os.path.exists(out_GC_content) or os.path.exists(out_GC_skew):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # convert
    fna2circos(args.input_file, out_GC_content, out_GC_skew, args.window_size)


if __name__ == "__main__":
    main()
