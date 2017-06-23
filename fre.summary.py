#!/usr/bin/env python
########################################################################
# Heteroplasmy liver lobes project:
# Summary of the heteroplasmy status of all analysed samples
#
# Copyright (c) 2017 Alexander Huebner <alexhbnr@gmail.com>
########################################################################

from __future__ import print_function
import argparse
import os
import re
import sys
import numpy as np
import pandas as pd


def main():
    ''' Summarises the heteroplasmy status of all analysed samples 
    '''
    outfile = open(Args['output'], "w")

    # Process heteroplasmy list
    hetlist = pd.read_csv(Args['input'], sep="\t", usecols=[0, 1, 10],
                          header=None, names=['sample', 'pos', 'maf'],
                          dtype={'sample': str, 'pos': int, 'maf': np.float64},
                          comment="#")
    hetlist['sample'] = hetlist['sample'].str.replace("_[a-z]+.log", "")
    hetlist['tissue'] = hetlist['sample'].str.replace("_[0-9]{3}", "")
    hetlist['id'] = hetlist['sample'].str.replace("^[A-Z][a-z]_|_[0-9]$", "")

    # Identify all IDs to be processed
    ids = sorted(list(set([re.sub("^[A-Z][a-z]_|.ssp$", "", f) if f.startswith("Bl")
                           else re.sub("^[A-Z][a-z]_|_[12].ssp$", "", f)
                           for f in os.listdir(path=Args['ssp'])])))
    # Iterate over all IDs
    for i in tqdm.tqdm(ids):
        tissues = [pd.read_csv(Args['ssp'] + "/Bl_{}.ssp".format(i), sep="\t",
                               usecols=[0, 2, 3, 4, 5, 6, 7, 8], index_col=[0]),
                   pd.read_csv(Args['ssp'] + "/Li_{}_1.ssp".format(i), sep="\t",
                               usecols=[0, 2, 3, 4, 5, 6, 7, 8], index_col=[0]),
                   pd.read_csv(Args['ssp'] + "/Li_{}_2.ssp".format(i), sep="\t",
                               usecols=[0, 2, 3, 4, 5, 6, 7, 8], index_col=[0])]
        for t in tissues:
            t['maf'] = t['second_count'] / t['cov']
        tissues[0].head()
        # Iterate overall positions of the MT genome
        for j in range(1, 16570):
            scenario = "".join([str(len(hetlist[(hetlist['id'] == i) & (hetlist['tissue'] == 'Bl') & (hetlist['pos'] == j)])),
                                str(len(hetlist[(hetlist['id'] == i) & (hetlist['tissue'] == 'Li_1') & (hetlist['pos'] == j)])),
                                str(len(hetlist[(hetlist['id'] == i) & (hetlist['tissue'] == 'Li_2') & (hetlist['pos'] == j)]))])
            first_allele = clean_alleles([t.at[j, 'first_allele'] for t in tissues],
                                         [t.at[j, 'cov'] for t in tissues],
                                         [t.at[j, 'first_count'] for t in tissues])
            second_allele = clean_alleles([t.at[j, 'second_allele'] for t in tissues],
                                          [t.at[j, 'cov'] for t in tissues],
                                          [t.at[j, 'second_count'] for t in tissues])
            other_allele = clean_alleles([t.at[j, 'other_allele'] for t in tissues],
                                         [t.at[j, 'cov'] for t in tissues],
                                         [t.at[j, 'other_count'] for t in tissues])
            mafs = [tissues[k].at[j, 'maf'] for k in range(0, 3)]
            mafs = [m if not np.isnan(m) else 0 for m in mafs]
            outfile.write("\t".join([str(i), str(j), scenario, first_allele, second_allele, other_allele] +
                                    ["{:.4f}".format(m) for m in mafs]) + "\n")

    outfile.flush()
    outfile.close()


def clean_alleles(alleles, coverages, counts):
    ''' Cleans alleles by removing duplicated bases and NA values
    '''
    allelefreqs = np.asarray(counts) / np.asarray(coverages)
    if allelefreqs[~np.isnan(allelefreqs)].shape[0] == 0:
        return ""
    else:
        max_allelefreq = np.nanmax(allelefreqs)
        maxid = np.argwhere(allelefreqs == max_allelefreq)[0][0]
        return alleles[maxid]


def status(*objs):
    ''' Allows to write status of programme to stderr
    '''
    print(*objs, file=sys.stderr)


Parser = argparse.ArgumentParser(description='Produces a FRE summary file')
Parser.add_argument('-i', '--input', required=True,
                    help='list of heteroplasmies')
Parser.add_argument('-s', '--ssp', required=True,
                    help='folder containing the SSP files per file')
Parser.add_argument('-o', '--output', required=True,
                    help='output FRE report')
Args = vars(Parser.parse_args())

if __name__ == '__main__':
    main()
