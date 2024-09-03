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


class BaseGffRecord(object):
    """ This is the base class of GFF record, all other GFF record types should
        inherite from this base GFF Record
    """
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attribute):
        """ Initialize a BaseGffRecord object, use standard GFF format specification.

        """
        self.seqname = seqname
        self.feature = feature
        self.source = source
        self.start = int(start)
        self.end = int(end)
        self.score = str(score)
        self.strand = strand
        self.frame = str(frame)
        self.attribute = attribute
        self.parentAttributeDict = None
        self.attribute_dict = self._get_attribute_dict()

    def _get_attribute_dict(self):
        attribute_dict = {}
        attri_list = self.attribute.strip().split(";")
        # need to tackle with "note=codon recognized: UUG; tRNA-Leu (CAA);"
        last_k, last_v = None, None
        for attr in attri_list:
            split_list = attr.strip().split("=")
            if len(split_list) == 2:
                if split_list[1] == "":
                    attribute_dict.update({split_list[0]:"None"})
                else:
                    attribute_dict.update({split_list[0]:split_list[1]})

                last_k, last_v = split_list
            else:
                assert len(split_list) == 1, "More than 2 parts were found !"
                if last_k:
                    attribute_dict[last_k] += ","+split_list[0]
                else:
                    raise Exception("No Key at all !")
        return attribute_dict

    def get_subattribute(self, subattribute):
        subattr = self.attribute_dict.get(subattribute, None)
        return subattr

    def get_parentSubAttribute(self, subattribute):
        subattr = self.parentAttributeDict.get(subattribute, None)
        return subattr

    def set_subattribute(self, subattribute, value):
        if self.attribute_dict.has_key(subattribute):
            self.attribute_dict[subattribute] = value
        else:
            self.attribute_dict.update({subattribute:value})
        # once set new value, should change self.attribute and self.attribute_dict
        self.attribute = ""
        for k, v in self.attribute_dict.iteritems():
            self.attribute += k+"="+v+";"
        self.attribute = self.attribute.rstrip(";")
        # update self.attribute_dict
        self.attribute_dict = self._get_attribute_dict()

    def __str__(self):
        return self.seqname+"\t"+self.source+"\t"+self.feature+"\t"+\
               str(self.start)+"\t"+str(self.end)+"\t"+self.score+"\t"+\
               self.strand+"\t"+self.frame+"\t"+self.attribute+"\n"


class BaseGffRecordParser(object):
    """ This class used to parse gff3 file, to generate BaseGffRecord instances
    """

    def __init__(self, handle_or_fileStr):
        self.handle = handle_or_fileStr
        self.genome_info = None

    def _line_parser(self):
        # judge file opened or not
        if not hasattr(self.handle, "read"):
            handle = open(self.handle, "r")
        else:
            handle = self.handle

        while True:
            line = handle.readline()

            if not line:
                # close file handle
                try:
                    handle.close()
                except Exception as e:
                    print(e)
                break

            else:
                if line.startswith("#"):
                    continue
                else:
                    gff_line_list = line.strip().split("\t")

                    gff_record = BaseGffRecord(*gff_line_list)
                    yield gff_record

    def __iter__(self):
        return self._line_parser()


def gff2circos(input_gff, out_plus_gene, out_minus_gene):
    """ given an input gff file, slice genes on the plus and minus strands,
        write out to out_plus_gene and out_minus_gene in bed format
    """
    with open(out_plus_gene, "w") as oh1, open(out_minus_gene, "w") as oh2:
        oh1.write("#chromosome\tstart\tend\tgene_name\n")
        oh2.write("#chromosome\tstart\tend\tgene_name\n")

        records = BaseGffRecordParser(input_gff)

        for rec in records:
            if rec.feature.lower() == "gene":
                chromosome = rec.seqname
                start = rec.start
                end = rec.end
                gene_name = rec.get_subattribute("Name")
                strand = rec.strand
                if strand == "+":
                    oh1.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+gene_name+"\n")
                else:
                    oh2.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+gene_name+"\n")


def main():


    # parser
    parser = argparse.ArgumentParser(description="create plus and minus strand\
                         gene content data for circos visualization")
    parser.add_argument("input_file", help="input gff file")
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
    out_plus_gene = os.path.join(args.out_folder, prefix+"_plus_gene.txt")
    out_minus_gene = os.path.join(args.out_folder, prefix+"_minus_gene.txt")

    if os.path.exists(out_plus_gene) or os.path.exists(out_minus_gene):
        if args.force:
            print("Warning: output file exists, will be overwriten!")
        else:
            print("Error: output file detected, please backup it at first")
            sys.exit(0)

    # convert
    gff2circos(args.input_file, out_plus_gene, out_minus_gene)


if __name__ == "__main__":
    main()
