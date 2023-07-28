#!/usr/bin/env python

# Copyright (C) 2020  Shengwei Hou
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


import sys, os
import click
import logging
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import OrderedDict
from collections import Counter



# logging
logfmt = "[%(asctime)s] [%(levelname)s] [%(name)s] %(message)s"
field_styles= {'asctime': {'color': 'green'}, 'levelname': {'bold': True, 'color': 'cyan'}, 
               'name': {'color': 'blue'}, 'username': {'color': 'magenta'}, 'hostname': {'color': 'yellow'}}
datefmt = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(stream=sys.stdout, format=logfmt, datefmt=datefmt)
logger = logging.getLogger('main')


def set_loglevel(loglevel):
    """
    :param loglevel: loglevel (str): minimum loglevel for emitting messages
    :return:
    """
    _level = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG}.get(loglevel, logging.DEBUG)

    logging.getLogger('main').setLevel(_level)





def calculate_aa_stats(input_faa, output_file, grep_str=None):

    aa_pct_dict = {}

    for rec in SeqIO.parse(input_faa, "fasta"):
        prot_name = rec.name
        prot_desc = rec.description
        prot_seq = str(rec.seq)
        analyzed_seq = ProteinAnalysis(prot_seq)
        #aa_count = analyzed_seq.count_amino_acids()
        aa_pct = analyzed_seq.get_amino_acids_percent()
        #print(aa_pct)
        if grep_str is not None:
            if grep_str not in prot_desc:
                continue
        aa_pct_dict[prot_name] = aa_pct

    aa_pct_df = pd.DataFrame.from_dict(aa_pct_dict, orient='index')
    aa_pct_df.to_csv(output_file, sep="\t", header=True, index=True, float_format="%.4f")
    print(aa_pct_df.mean(axis=0))


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('input_faa')
@click.option('-g', '--grep_str',
              type=click.STRING,
              default=None,
              show_default=True,
              help="string used to grep sequence, such as using 'ribosomal protein' to select only ribosomal proteins")
@click.option('-o', '--output_dir',
              type=click.STRING,
              default='./',
              show_default=True,
              help="output directory")
@click.option('-p', '--prefix',
              type=click.STRING,
              default=None,
              show_default=True,
              help="output prefix")
@click.option('-f', '--force',
              is_flag=True,
              help="force to run")
@click.option('-l', '--loglevel', 
              default='info', 
              show_default=True,
              type=click.Choice(['critical', 'error', 'warning', 'info', 'debug']), 
              help="logging level")
def main(input_faa, grep_str, output_dir, prefix, force, loglevel):

    # logging
    set_loglevel(loglevel)

    # output handeling
    if not prefix:
        prefix = os.path.splitext(os.path.split(input_faa)[-1])[0]
    output_prefix = os.path.join(output_dir, prefix)
    output_file = output_prefix + "_aa_stats.tsv"
    if os.path.exists(output_file):
        if not force:
            logger.error("{f} exists, use --force to overwrite".format(f=output_file))
            sys.exit(0)

    # calculate aa stats
    calculate_aa_stats(input_faa, output_file, grep_str)


if __name__=='__main__':
    main()
