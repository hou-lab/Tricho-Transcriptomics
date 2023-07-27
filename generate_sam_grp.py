#!/usr/bin/env python


# <generate_sam_grp.py, do segemehl.x alignment then produce grp file.>
# Copyright (C) <2016>  <Shengwei Hou> <housw2010'at'gmail'dot'com>
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


import sys
import os
import re
import subprocess
from Bio import SeqIO
import numpy as np
import argparse


class SamRecord(object):
    """ class to represent the SAM record in segemehl sam file
    """

    __slots__ = ["qname","shortname", "oriStrand","weight", "mapStrand",
                 "rname", "startPos", "mapq", "cigar",
                 "rnext", "pnext", "tlen", "seq","qual", "option"]

    def __init__(self, qname, shortname, oriStrand, weight, mapStrand, rname,
                 startPos, mapq, cigar, rnext, pnext, tlen,
                 seq, qual, option):

        self.qname = qname
        self.shortname = shortname
        self.oriStrand = oriStrand
        self.weight = weight
        self.mapStrand = mapStrand
        self.rname = rname
        self.startPos = startPos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.option = option
        self._update_weight()

    def _update_weight(self):
        """ divide the weight according to # of hit positions
        """
        _repeat_time = int(self.option['NH'])
        self.weight /= float(_repeat_time)

    def __str__(self):
        return self.qname + " mapped to "+self.mapStrand + \
        " strand, locates in <" + self.rname + "> starts at <" + \
        str(self.startPos)+">\n"



class SamParser(object):
    """ class to parse the segemehl SAM record from a open file handle
    """

    def __init__(self, handle_or_fileStr):

        self.handle = handle_or_fileStr

    def _lines(self):

        # judge file opened or not
        if not hasattr(self.handle, "read"):
            handle = open(self.handle, "r")
        else:
            handle = self.handle

        while True:
            line = handle.readline()
            if not line:
                handle.close()
                break
            else:
                if line.startswith("@"):
                    continue
                else:
                    line = line.strip()
                    yield line

    def _parse_SAM(self, inline):
        for line in inline:
            sam_line = line.strip().split("\t")
            fullname = sam_line[0].replace("_", " ")
            qname = " ".join(item for item in fullname.split(" ")[0:2])
            shortname = fullname.split(" ")[0]
            try:
                oriStrand = int(fullname.split(" ")[1][0])
            except Exception as e:
                oriStrand = 1
            if "weight|" in fullname:
                weight = float(fullname.strip().split("|")[1])
            else:
                weight = 1

            # the strand info that reads mapped in reference genome
            if int(sam_line[1]) == 16 or int(sam_line[1]) == 272:
                mapStrand = "-"
            elif int(sam_line[1]) == 0 or int(sam_line[1]) == 256:
                mapStrand = "+"
            else:
                continue 

            # cigar info and use it to get map info
            cigar = sam_line[5]
            # use left most coordinate system
            startPos = int(sam_line[3])
            rname = sam_line[2]
            mapq = int(sam_line[4])
            rnext = sam_line[6]
            pnext = int(sam_line[7])
            tlen = int(sam_line[8])
            seq = sam_line[9]
            qual = sam_line[10]

            # put all additional fields into option
            option = {}

            for item in sam_line[11:]:
                parts = item.strip().split("\t")
                for part in parts:
                        new_part = part.strip().split(":")
                        option.update({new_part[0]: new_part[2]})

            sam_record = SamRecord(qname, shortname, oriStrand, weight, mapStrand,
                                   rname, startPos, mapq, cigar,
                                   rnext, pnext, tlen, seq, qual, option)

            yield sam_record

    def __iter__(self):
        return self._parse_SAM(self._lines())


def do_segemehl_align(input_fasta, input_genome, out_folder, identity=75):
    """ function to do segemehl alignment
    """

    input_basename = os.path.splitext(os.path.split(input_fasta)[-1])[0]
    genome_basename = os.path.splitext(os.path.split(input_genome)[-1])[0]

    # mkdir if out_folder not exits
    if not os.path.exists(out_folder):
        subprocess.check_call(["mkdir", "-p", out_folder])

    out_file = os.path.join(out_folder, input_basename+"_mapped_to_"+genome_basename+"_id"+str(identity)+".sam")
    out_notmatched = os.path.join(out_folder, input_basename+"_mapped_to_"+genome_basename+"_id"+str(identity)+"_non_matched.txt")


    try:
        subprocess.check_call(["segemehl.x",
                               "--database", input_genome,
                               "-x", genome_basename+"index",
                               "--query", input_fasta,
                               "--outfile", out_file,
                               "--minsize", "20",
                               "--threads", "20",
                               "-u", out_notmatched,
                               "--extensionpenalty", "4",
                               "--dropoff", "8",
                               "--accuracy", str(identity),
                               "--minsplicecover", "80",
                               "--minfraglen", "20",
                               "--hitstrategy", "1"
                              ])
    except Exception as e:
        print(e)

    return out_file, out_notmatched


def createFolder(path2fileName):
    """ function to create folder under given path
    """
    fileExist = False
    path2dir = path2fileName
    try:
        os.mkdir(path2dir)
    except OSError as e:
        # File already exists
        print("\nWarning: The original '%s' folder already exists, will use it !"%path2dir)
        fileExist = True
    return fileExist


def cigar_parser(cigar, retStr=False):
    """
        Op    BAM    Description
    ---------------------------------------------
        M :    0     alignment Match (can be sequence match or mismatch)
        I :    1     Insert to the reference
        D :    2     Deletion from the refrence
        N :    3     skipped region from the reference
        S :    4     Soft clipping in Field SEQ (sequence present)
        H :    5     Hard clipping in Field SEQ (sequence NOT present)
        P :    6     Padding (silent deletionfrom padded reference)
        = :    7     sequence MATCH
        X :    8     sequence MISMATCH

        parse cigar string, give operator and operation times subsequently

        >>> cigar_parser("3S17M8D4M9I3H")

        [('S', '3'), ('M', '17'), ('D', '8'), ('M', '4'), ('I', '9'), ('H', '3')]

        >>> cigar_parser("3S17M8D4M9I3H", retStr=True)
        "SSSMMMMMMMMMMMMMMMMMDDDDDDDDMMMMIIIIIIIIIHHH"

    """

    result = []
    n = ''
    for c in cigar:
        # concatenate digits
        if c.isdigit():
            n += c
        # count operator and operation times
        elif c in ["M", "I", "D", "N"]:
            if n == "":
                raise ValueError("End of CIGAR string reached, but an operator was expected")
            result.append((c, n))
            n = ''
        # negelect other operators not exist in above list
        else:
            n = ''
            continue

    if not retStr:
        return result
    else:
        ret = ""
        for tup in result:
            ret += tup[0]*int(tup[1])
        return ret


def generate_grp_file(genome_length, sam_file, results_folder):
    """ function to generate grp file from input sam file
    """
    basename = os.path.splitext(os.path.split(sam_file)[-1])[0]
    grp_file = os.path.join(results_folder, basename+".grp")

    # initialize array as long as genome length
    pos_array = np.arange(genome_length)
    fwd_cov = np.zeros(genome_length)
    rev_cov = np.zeros(genome_length)
    fwd_tss = np.zeros(genome_length)
    rev_tss = np.zeros(genome_length)

    # parse sam file
    samFile = SamParser(sam_file)

    # strand info
    strand = None

    for record in samFile:
        full_cigar = cigar_parser(record.cigar, retStr=True)
        # fwd reads, generate tss file
        if int(record.oriStrand) == 1:

            # make sure all reads are from FV library
            if not strand:
                strand = 1
            else:
                assert strand == 1, "mixed reads origin in this sam file !!"

            # fwd reads mapped to + strand
            if record.mapStrand == "+":

                # update tss
                fwd_tss[record.startPos-1] += record.weight

                # update cov
                i = record.startPos -2
                for char in full_cigar:
                    # for match or deletion in reference, count coverage
                    if char == "M" or char == "D":
                        i += 1
                        fwd_cov[i] += record.weight
                    # for insert in reference, do nothing
                    elif char == "I":
                        continue
                    # for intron in reference, walk through
                    elif char == "N":
                        i +=1
                    else:
                        print("unexpected char in cigar, ", char)
                        continue

            # fwd reads mapped to - strand
            else:
                assert record.mapStrand == "-", "mapStrand can only be + or - !"

                # update cov
                i = record.startPos -2
                for char in full_cigar:
                    # for match or deletion in reference, count coverage
                    if char == "M" or char == "D":
                        i += 1
                        rev_cov[i] += record.weight
                    # for insert in reference, do nothing
                    elif char == "I":
                        continue
                    # for intron in reference, walk through
                    elif char == "N":
                        i += 1
                    else:
                        print("unexpected char in cigar, ", char)
                        continue

                # update tss
                rev_tss[i] += record.weight

        # rev reads, now also generate tss file, but will colored into gray finally
        else:
            assert int(record.oriStrand) ==2, "oriStrand can only be 1 or 2 !"

            # make sure all reads are from RV library
            if not strand:
                strand = 2
            else:
                assert strand == 2, "mixed reads origin in this sam file !!"

            # rev reads mapped to + strand, comes from reverse gene
            if record.mapStrand == "+":

                # update cov
                i = record.startPos -2
                for char in full_cigar:
                    # for match or deletion in reference, count coverage
                    if char == "M" or char == "D":
                        i += 1
                        rev_cov[i] += record.weight
                    # for insert in reference, do nothing
                    elif char == "I":
                        continue
                    # for intron in reference, walk through
                    elif char == "N":
                        i += 1
                    else:
                        print("unexpected char in cigar, ", char)
                        continue

                # update tss
                rev_tss[i] += record.weight

            # rev reads mapped to - strand, comes from forward gene
            else:
                assert record.mapStrand == "-", "mapStrand can only be + or - !"

                # update tss
                fwd_tss[record.startPos-1] += record.weight

                # update cov
                i = record.startPos -2
                for char in full_cigar:
                    # for match or deletion in reference, count coverage
                    if char == "M" or char == "D":
                        i += 1
                        fwd_cov[i] += record.weight
                    # for insert in reference, do nothing
                    elif char == "I":
                        continue
                    # for intron in reference, walk through
                    elif char == "N":
                        i +=1
                    else:
                        print("unexpected char in cigar, ", char)
                        continue

    pos_array += 1

    # make rev to minus
    rev_cov *= -1
    rev_tss *= -1

    # fwd reads have cov and tss info, tss has colors
    if strand == 1:
        with open(grp_file, "w") as oh:
            # to customize the line color, we must use space delimited format
            oh.write("# BASE fwd_coverage fwd_tss rev_coverage rev_tss\n")
            oh.write("# colour 0:204:0 255:0:0 0:0:255 255:0:0\n")
            for pos, fwd_c, fwd_t, rev_c, rev_t in zip(pos_array, fwd_cov, fwd_tss, rev_cov, rev_tss):
                oh.write(str(pos)+" %.2f %.2f %.2f %.2f\n"%(fwd_c, fwd_t, rev_c, rev_t))

    # rev reads have also cov and tss info, but this tss are not so useful,
    # will be colored as gray
    else:
        assert strand ==2, "strand can only be 1 or 2 !"
        with open(grp_file, "w") as oh:
            # to customize the line color, we must use space delimited format
            oh.write("# BASE fwd_coverage fwd_tss rev_coverage rev_tss\n")
            oh.write("# colour 0:204:0 160:160:160 0:0:255 160:160:160\n")
            for pos, fwd_c, fwd_t, rev_c, rev_t in zip(pos_array, fwd_cov, fwd_tss, rev_cov, rev_tss):
                oh.write(str(pos)+" %.2f %.2f %.2f %.2f\n"%(fwd_c, fwd_t, rev_c, rev_t))


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="do segemehl.x alignment then produce grp file")
    parser.add_argument("in_reads", help="input reads in fasta format")
    parser.add_argument("in_ref", help="input genome reference in fasta format")
    parser.add_argument("identity", help="minimum identity for segemehl alignment")
    parser.add_argument("-o", "--out_folder", default="sam_grp", help="output folder to put results")
    args = parser.parse_args()

    input_fasta = args.in_reads
    input_genome = args.in_ref
    identity = args.identity
    results_folder = args.out_folder

    # ---------------------------------

    # get genome length
    genomeFastaRecord = SeqIO.read(input_genome, "fasta")
    genome_length = len(genomeFastaRecord.seq)

    # ----------------------------------

    # create a folder to put results generated by segemehl
    segemehl_results_folder = os.path.join(results_folder, "after_segemehl")
    createFolder(segemehl_results_folder)

    # align input_fasta to input_genome using segemehl
    sam_file, not_matched_file = do_segemehl_align(input_fasta, input_genome,
                                                   segemehl_results_folder, identity)

    # ----------------------------------

    # create a folder to put grp results
    grp_results_folder = os.path.join(results_folder, "grp_files")
    createFolder(grp_results_folder)

    # generate grp file
    generate_grp_file(genome_length, sam_file, grp_results_folder)


if __name__ == "__main__":
    main()
