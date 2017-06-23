#!/usr/bin/env python
###############################################################################
# Heteroplasmy liver lobes project:
# Filters all heteroplasmies that might have suffered a cross-contamination
# during processing in the wet lab
#
# Copyright (c) 2017 Alexander Huebner <alexhbnr@gmail.com>
###############################################################################

from __future__ import print_function
import argparse
from itertools import cycle
import pandas as pd


def main():
    ''' Filters the list of heteroplasmies for the presence of cross-contaminants
    '''
    # Read heteroplasmy list
    hetlist = pd.read_csv(Args['input'], sep="\t", header=None, comment="#")
    hetlist['id'] = hetlist[0].str.replace("_[a-z]+.log", "")
    hetlist = hetlist.set_index(['id'])

    # Read list of contaminants
    for contfn in Args['cont']:
        contlist = pd.read_csv(contfn, sep="\t", header=None, usecols=[0, 1])
        id_list = list(set(contlist[0]))
        tissue_list = list(set(contlist[1]))
        contaminants = ["_".join([t, str(i)]) for i, t in zip(id_list, cycle(tissue_list))]
        hetlist.drop(contaminants, inplace=True)

    # Write updated list to file
    hetlist.to_csv(Args['output'], sep="\t", header=False, index=False, float_format="%.3f")


Parser = argparse.ArgumentParser(description='Filter the list of heteroplasmies' +
                                 ' for the presence of cross-contaminated samples')
Parser.add_argument('-i', '--input', required=True,
                    help='list of heteroplasmies')
Parser.add_argument('-c', '--cont', required=True, nargs="+",
                    help='list of cross-contaminants')
Parser.add_argument('-o', '--output', required=True,
                    help='filtered heteroplasmy list')
Args = vars(Parser.parse_args())


if __name__ == '__main__':
    main()
