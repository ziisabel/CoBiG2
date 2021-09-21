#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import mainFunctions as mf


'''Generates text files (.tab) to use as input in MapMan to perform pathway analysis'''


def writeTairsFC(allDegs, degsTair, DE_dir, outdir):

    """
    Generates .tab files with TAIRs homologs to DEGs and respective log2FC to use with MapMan

    Parameters:
        list of all DEGs detected through DE analysis
        dictionary with genes as keys and the respective mapped TAIRs (blast) as values
        directory with all results files from DE analysis (.csv)
    
    Output: .tab files with TAIRs and log2FC, ordered by log2FC (most to least up- and down-regulated)
    """

    # saves DE results into a dictionary with comparisons as keys and log2FC as values
    for file in [f for f in os.listdir(DE_dir) if not f.startswith('.')]:
        with open(DE_dir + file, 'r') as infile:
            infile.readline()
            DE_results = {line.split(',')[0].replace('\"',''): float(line.split(',')[2]) for line in infile.readlines() if line.split(',')[0].replace('\"','') in allDegs}
        infile.close()

        # creates a dictionary grouped by regulation (up- and down-regulated DEGs) with DE results ordered by log2FC
        degsFC = {}
        degsFC['up'] = {k: v for k, v in sorted(DE_results.items(), key=lambda item: item[1], reverse=True) if v > 0}
        degsFC['down'] = {k: v for k, v in sorted(DE_results.items(), key=lambda item: item[1], reverse=True) if v < 0}

        # writes 2 .tab file per comparison (up- and down-regulated DEGs) with the TAIRs and log2FC respective to its DEGs
        for regulation in degsFC:
            file_out = outdir + file.strip('.csv') + '_' + regulation + '.txt'
            with open(file_out, 'w') as outfile:
                for gene in degsFC[regulation]:
                    if gene in degsTair:
                        outfile.write(degsTair[gene] +'\t'+ degsFC[regulation][gene] +'\n')
            outfile.close()


writeTairsFC(mf.getAllDegsSet(sys.argv[1]), mf.getDegsTair(sys.argv[2]), sys.argv[3], sys.argv[4])