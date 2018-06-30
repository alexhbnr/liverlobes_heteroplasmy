#!/usr/bin/env python
########################################################################
# Heteroplasmy liver lobes project:
# Summary of the heteroplasmy status of all analysed samples
#
# Copyright (c) 2017 Alexander Huebner <alexhbnr@gmail.com>
########################################################################

from __future__ import print_function
import argparse
import tqdm
import numpy as np
import pandas as pd


def main():
    ''' Summarises the heteroplasmy status of all analysed samples
    '''
    # Overwrite previous file content due to append mode later
    outfile = open(Args['output'], "w")
    outfile.close()

    # Process heteroplasmy list
    hetlist = pd.read_csv(Args['input'], sep="\t", usecols=[0, 1, 10],
                          header=0, names=['sample', 'pos', 'maf'],
                          dtype={'sample': str, 'pos': int, 'maf': np.float64},
                          comment="#")
    hetlist['tissue'] = hetlist['sample'].str.extract("S[0-9]+_([A-Za-z]+_*[0-9]*)", expand=True)
    hetlist['id'] = hetlist['sample'].str.extract("S([0-9]+)_[A-Za-z]+_*[0-9]*", expand=True)

    # Identify all IDs to be processed
    ids = hetlist['id'].unique()
    tissue_cats = ["Bl", "Li_1", "Li_2"]

    # Iterate over all IDs
    for i in tqdm.tqdm(ids):
        tissues = [pd.read_csv(Args['ssp'] + "/S{}_Bl.ssp".format(i), sep="\t",
                               usecols=[0, 2, 3, 4, 5, 6, 7, 8], index_col=[0]),
                   pd.read_csv(Args['ssp'] + "/S{}_Li_1.ssp".format(i), sep="\t",
                               usecols=[0, 2, 3, 4, 5, 6, 7, 8], index_col=[0]),
                   pd.read_csv(Args['ssp'] + "/S{}_Li_2.ssp".format(i), sep="\t",
                               usecols=[0, 2, 3, 4, 5, 6, 7, 8], index_col=[0])]
        # Clean up NaN and add additinal columns to SSP files
        for j, t in enumerate(tissues):
            t.loc[:, ("first_allele", "second_allele", "other_allele")] = t.loc[:, ("first_allele", "second_allele", "other_allele")].fillna("", axis=1)
            t.loc[:, ("first_count", "second_count", "other_count")] = t.loc[:, ("first_count", "second_count", "other_count")].fillna(0, axis=1)
            t['maf'] = t['second_count'] / t['cov']
            t['het'] = 0
            t.loc[np.isin(t.index.values, hetlist.loc[(hetlist['id'] == i) & (hetlist['tissue'] == tissue_cats[j])]["pos"].tolist()), 'het'] = 1

        # Construct FRE file
        # Heteroplasmy scenario
        heteroplasmy = pd.DataFrame(dict(Bl=tissues[0]['het'],
                                         Li_1=tissues[1]['het'],
                                         Li_2=tissues[2]['het']))
        heteroplasmy['heteroplasmy'] = heteroplasmy.apply(sum, axis=1) > 0
        heteroplasmy['scenario'] = heteroplasmy.apply(lambda row: "".join([str(r) for r in row[:3]]), axis=1)
        if heteroplasmy[heteroplasmy['heteroplasmy']].shape[0] > 0:
            # Calculate MAF
            maf = pd.DataFrame(dict(Bl=tissues[0]['maf'],
                                    Li_1=tissues[1]['maf'],
                                    Li_2=tissues[2]['maf'])) \
                    .loc[heteroplasmy['heteroplasmy']]
            # Check the major alleles: if different between liver and blood, use blood as base line
            allele_discrepancy_mask = pd.DataFrame(dict(major_Bl=tissues[0]['first_allele'],
                                                        major_Li_1=tissues[1]['first_allele'],
                                                        major_Li_2=tissues[2]['first_allele'],
                                                        minor_Bl=tissues[0]['second_allele'],
                                                        minor_Li_1=tissues[1]['second_allele'],
                                                        minor_Li_2=tissues[2]['second_allele'])) \
                                        .loc[heteroplasmy['heteroplasmy']]
            allele_discrepancy_mask['nMajor'] = allele_discrepancy_mask.apply(lambda row: row[:3].nunique(), axis=1)
            allele_discrepancy_mask['nMinor'] = allele_discrepancy_mask.apply(lambda row: row[3:6].nunique(), axis=1)
            allele_discrepancy_mask['FlipMajor'] = (allele_discrepancy_mask['nMajor'] > 1) & (allele_discrepancy_mask['nMinor'] > 1)
            allele_discrepancy_mask['minorAllele'] = allele_discrepancy_mask.loc[:, ('major_Bl', 'major_Li_1', 'major_Li_2',
                                                                                     'minor_Bl', 'minor_Li_1', 'minor_Li_2', 'FlipMajor')] \
                                                        .join(heteroplasmy.loc[heteroplasmy['heteroplasmy'], ('Bl', 'Li_1', 'Li_2')]) \
                                                        .apply(minorAllele, axis=1)
            flip_sites = allele_discrepancy_mask[allele_discrepancy_mask['FlipMajor']].index.values
            if (allele_discrepancy_mask.loc[flip_sites, 'major_Bl'] !=
                    allele_discrepancy_mask.loc[flip_sites, 'major_Li_1']).any():
                maf.loc[flip_sites, "Li_1"] = (tissues[1].loc[flip_sites, "first_count"] /
                                               tissues[1].loc[flip_sites, "cov"])
            if (allele_discrepancy_mask.loc[flip_sites, 'major_Bl'] !=
                    allele_discrepancy_mask.loc[flip_sites, 'major_Li_2']).any():
                maf.loc[flip_sites, "Li_2"] = (tissues[2].loc[flip_sites, "first_count"] /
                                               tissues[2].loc[flip_sites, "cov"])

            # Prepare output
            fre = maf.join(heteroplasmy['scenario'])
            fre['id'] = i
            fre['majorAllele'] = tissues[0].loc[heteroplasmy['heteroplasmy'], "first_allele"]
            fre['minorAllele'] = allele_discrepancy_mask['minorAllele']
            fre = fre.reset_index()[['id', 'pos', 'scenario', 'majorAllele', 'minorAllele', 'Bl', 'Li_1', 'Li_2']]
            fre.to_csv(Args['output'], mode="a", sep="\t", header=False, float_format="%.3f", index=False)


def minorAllele(site):
    ''' Determine the correct minor allele for the FRE file. If blood has a
    minor allele, then its minor allele is taken. If blood doesn't have one, we
    will check, which liver sample is heteroplasmic for the site and use their
    minor allele instead.
    '''
    if site[3] not in ["", "I", "D"]:
        minor_allele = site[3]
    else:
        if site[6]:  # flip allele
            het_minor_alleles = list(set([site[i] for i in range(1, 3) if site[i+7] == 1 and site[i] != site[0]]))
            if len(het_minor_alleles) == 1:
                minor_allele = het_minor_alleles[0]
            else:
                minor_allele = "?"
        else:
            het_minor_alleles = list(set([site[i] for i in range(4, 6) if site[i+4] == 1 and site[i] != site[0]]))
            if len(het_minor_alleles) == 1:
                minor_allele = het_minor_alleles[0]
            else:
                minor_allele = "?"
    return minor_allele


if __name__ == '__main__':
    Parser = argparse.ArgumentParser(description='Produces a FRE summary file')
    Parser.add_argument('-i', '--input', required=True,
                        help='list of heteroplasmies')
    Parser.add_argument('-s', '--ssp', required=True,
                        help='folder containing the SSP files per file')
    Parser.add_argument('-o', '--output', required=True,
                        help='output FRE report')
    Args = vars(Parser.parse_args())

    main()
