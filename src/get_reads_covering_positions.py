#!/usr/bin/env python
"""Get all reads covering positions file"""
########################################################################
# File: get_reads_covering_positions.py
#  executable: get_reads_covering_positions.py
#
# Author: Andrew Bailey
# History: 05/19/20 Created
########################################################################

import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from py3helpers.utils import time_it
from subprocess import check_call


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--bam', '-d', action='store',
                        dest='bam', required=True, type=str, default=None,
                        help="Bam file to look through")
    parser.add_argument('--positions_file', '-p', action='store',
                        dest='positions_file', required=True, type=str, default=None,
                        help="Path to base positions file")
    parser.add_argument('--output_file', '-o', action='store',
                        dest='output_file', required=True, type=str, default=None,
                        help="Path to new bam file")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    assert os.path.exists(args.bam), "{} does not exist".format(args.bam)
    assert os.path.exists(args.positions_file), "{} does not exist".format(args.positions_file)
    step_size = 10000
    step_number = 0
    with open(args.positions_file, "r") as fh:
        while True:
            execute = "samtools view -b -o {} {}".format(args.output_file+str(step_number), args.bam)
            positions = " "
            counter = 0
            for line in fh:
                split_line = line.split()
                chromosome = split_line[0]
                position = split_line[1]
                positions += chromosome+":"+position+"-"+position+" "
                counter += 1
                if counter > 10000:
                    break
            if positions == " " or step_number > 10:
                break
            execute += positions
            check_call(execute.split())
            step_number += 1


if __name__ == '__main__':
    print(time_it(main)[1])
