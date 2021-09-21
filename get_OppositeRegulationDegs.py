#!/usr/bin/env python3
# coding: utf-8

import mainFunctions as mf
import sys


'''Creates two .txt tab-delimited files (up and down) with lists of oppositely regulated DEGs shared between two comparisons.'''


def getAnnotations(filename):

    """
    Creates a dictionary with genes annotations, with gene/transcript IDs as keys and gene and protein names as values
    
    Parameter: tab-delimited filename containing genes annotations, with at least gene ID, gene name and protein name

    Return: dictionary with genes annotations, with gene/transcript IDs as keys and gene and protein names as values
    """

    annotations = {}

    for line in mf.readFile(filename, True):
        line = line.split('\t')
        annotations[line[0]] = {'gene': line[2], 'protein': line[3]}

    return annotations


def writeOppositeDegs(comparisonA, comparisonB, outdir):

    """
    Creates two .txt tab-delimited files (up and down) with lists of oppositely regulated DEGs shared between two comparisons
    
    Parameter:
        comparisons to compared
        output directory
    """

    allDegs = mf.getAllDegsDict('../DEGs/', '\t', 3)

    opposites = {}

    for comparison in allDegs:
        opposites[comparison] = {'up': {}, 'down': {}}

    for geneid in allDegs[comparisonA]:
        if geneid in allDegs[comparisonB]:
            if '-' in allDegs[comparisonA][geneid]:
                if '-' not in allDegs[comparisonB][geneid]:
                    opposites[comparisonA]['down'][geneid] = {'gene': annotations[geneid]['gene'], 'protein': annotations[geneid]['protein']}
            else:
                if '-' in allDegs[comparisonB][geneid]:
                    opposites[comparisonA]['up'][geneid] = {'gene': annotations[geneid]['gene'], 'protein': annotations[geneid]['protein']}

    for regulation in opposites[comparisonA]:
        with open(outdir + comparisonA + '_' + comparisonB + '_' + regulation + '.txt', 'w') as outfile:
            outfile.write('Trinity ID' +'\t'+ 'Gene' +'\t'+ 'Protein name' +'\n')
            for geneid in opposites[comparisonA][regulation]:
                outfile.write(geneid +'\n')
        outfile.close()


annotations = getAnnotations(sys.argv[1])
writeOppositeDegs(sys.argv[2], sys.argv[2], sys.argv[3])