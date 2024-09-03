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
import pandas as pd

def DE2circos(input_file, prefix,
              chromosome_col, start_col, end_col, strand_col, geneID_col,
              number_col_names):
    """ given an input annotatated DE file, slice columns to compose a circos
        bed file for fwd and rev strand, which contains

        chromosome start end  value

    """
    # slice values for each given names
    for value_name in number_col_names:
        with open(input_file, "r") as ih, open(prefix+"_plus_"+value_name+".txt", "w") as oh1, \
             open(prefix+"_minus_"+value_name+".txt", "w") as oh2:
            oh1.write("#chromosome\tstart\tend\t%s\n"%value_name)
            oh2.write("#chromosome\tstart\tend\t%s\n"%value_name)
            df = pd.read_table(input_file, header=0)
            for i in df.index:
                chromosome = df.loc[i, chromosome_col]
                start = df.loc[i, start_col]
                end = df.loc[i, end_col]
                strand = df.loc[i, strand_col]
                value = df.loc[i, value_name]
                if strand == "+":
                    oh1.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+str(value)+"\n")
                else:
                    oh2.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+str(value)+"\n")

    # slice up and down regulated genes
    with open(input_file, "r") as ih, open(prefix+"_DE_plus_up.txt", "w") as oh1, \
             open(prefix+"_DE_minus_up.txt", "w") as oh2, \
             open(prefix+"_DE_plus_down.txt", "w") as oh3, \
             open(prefix+"_DE_minus_down.txt", "w") as oh4:
            oh1.write("#chromosome\tstart\tend\tgeneName\n")
            oh2.write("#chromosome\tstart\tend\tgeneName\n")
            oh3.write("#chromosome\tstart\tend\tgeneName\n")
            oh4.write("#chromosome\tstart\tend\tgeneName\n")
            df = pd.read_table(input_file, header=0)
            for i in df.index:
                chromosome = df.loc[i, chromosome_col]
                start = df.loc[i, start_col]
                end = df.loc[i, end_col]
                strand = df.loc[i, strand_col]
                geneName = df.loc[i, geneID_col]
                logFC = df.loc[i, number_col_names[0]]
                padj = df.loc[i, number_col_names[1]]
                if float(logFC) > 0 and float(padj) <= 0.01:
                    if strand == "+":
                        oh1.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+geneName+"\n")
                    else:
                        oh2.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+geneName+"\n")
                elif float(logFC) < 0 and float(padj) <= 0.01:
                    if strand == "+":
                        oh3.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+geneName+"\n")
                    else:
                        oh4.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+geneName+"\n")
                else:
                    continue


def main():

    # parser
    parser = argparse.ArgumentParser(description="create plus and minus strand\
                         DE gene data for circos visualization")
    parser.add_argument("input_file", help="input annotatated DE table")
    parser.add_argument("-c", "--chromosome_col_name", required=False, default="seqname",
                              type=str,
                              help="column name where chromosome is located, default=seqname")
    parser.add_argument("-s", "--start_col_name", required=False, default="start", type=str,
                              help="column name where gene start is located, default=start")
    parser.add_argument("-e", "--end_col_name", required=False, default="end", type=str,
                              help="column name where gene end is located, default=end")
    parser.add_argument("-d", "--strand_col_name", required=False, default="strand", type=str,
                              help="column name where gene strand is located, default=strand")
    parser.add_argument("-g", "--geneID_col_name", required=False, default="geneID", type=str,
                              help="column name where geneID is located, default=geneID")
    parser.add_argument("-n", "--number_col_names", required=False, default=["log2FoldChange", "padj"], 
                              type=str, nargs="+", 
                              help="column names where gene associated numbers are located,"\
                              "like foldchange or padj, the first two names should be "\
                              "log2FoldChange and padj in order to slice DE genes. default=['log2FoldChange', 'padj']")
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
    prefix = os.path.join(args.out_folder, prefix)
    
    # convert
    DE2circos(args.input_file, prefix,
              args.chromosome_col_name, args.start_col_name, args.end_col_name,
              args.strand_col_name, args.geneID_col_name, args.number_col_names)


if __name__ == "__main__":
    main()
