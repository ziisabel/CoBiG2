#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import mainFunctions as mf


'''Creates multiple .txt tab-delimited files with lists of knockout (KO) genes in a sample, relative to its test or control.'''


def getGeneCounts(counts_file):

    """
    Creates a dictionary with gene IDs, the sample that expressed it and its gene counts

    Parameter: fullname of the file containing as samples' gene counts

    Return: a dictionary with gene IDs as keys and a dictionary as value with sample name as key and its gene count as value
    """

    with open(counts_file, 'r') as infile:
        header = infile.readline().split('\t')
        content = [i.strip('\n') for i in infile.readlines()]
    infile.close()

    geneCounts = {}

    for line in content:
        line = line.split('\t')
        geneid = line[0]
        geneCounts[geneid] = {}

        for attribute in header[1:]:
            # get clean sample name from filename
            sample = attribute.split('_')
            if len(sample) > 3:
                sample = sample[0] + '_' + sample[2]
            else:
                sample = sample[0] + '_' + sample[1]

            # sum the gene counts from all replicates of a sample
            if sample not in geneCounts[geneid]:
                geneCounts[geneid][sample] = int(float(line[header.index(attribute)]))
            else:
                geneCounts[geneid][sample] += int(float(line[header.index(attribute)]))
        
    return geneCounts


def getAllDegs(indir):

    """
    Reads DEGs from all DE results files and saves them as a dictionary with comparisons (DE) as keys and the respective list of DEGs as values

    Parameter:
        directory containing all DE results files (.csv)

    Return: dictionary with comparisons (DE) as keys and the respective list of DEGs as values
    """

    allDegs = {}

    for filename in [f for f in os.listdir(indir) if not f.startswith('.')]:

        degs = []

        for line in mf.readFile(indir + filename, True):
            line = line.strip('\n').split('\t')
            # add an entry for each gene to the degs dictionary {gene: log2FC}
            degs.append(line[0])
        
        # add the degs dictionary as value to the comparison dictionary {comparison: {degs}}
        allDegs[filename.strip('.tab')] = degs

    return allDegs


def getDegsCounts(geneCounts, allDegs):

    """"
    Creates a dictionary with DE comparison names as keys and the respective lists of KO gene IDs as values

    Parameter:
        output of the getGeneCounts() function
        output of the getAllDegs() function

    Return: dictionary with DE comparison names as keys and the respective lists of KO gene IDs as values
    """

    knockout_degs = {}

    for comparison in allDegs:
        knockout_degs[comparison] = set()

        test = comparison.split('_vs_')[0]
        control = comparison.split('_vs_')[1]
        
        for geneid in geneCounts:
            if geneid in allDegs[comparison]:
                
                if geneCounts[geneid][control] == 0:
                    if geneCounts[geneid][test] > 0:
                        knockout_degs[comparison].add(geneid)
                if geneCounts[geneid][test] == 0:
                    if geneCounts[geneid][control] > 0:
                        knockout_degs[comparison].add(geneid)
    
    return knockout_degs


def writeKnockoutDegs(knockout_degs, outdir):

    """
    Writes multiple .txt tab-delimited files (one per DE comparison) containing a list of the respective knockout genes (one gene per line)

    Parameters:
        output of the getDegsCounts() function
        fullname of the output directory
    """

    for comparison in knockout_degs:
        with open(outdir + comparison + '.txt', 'w') as outfile:
            for geneid in knockout_degs[comparison]:
                outfile.write(geneid + '\n')


knockout_degs = getDegsCounts(getGeneCounts(sys.argv[1]), getAllDegs(sys.argv[2]))
writeKnockoutDegs(knockout_degs, sys.argv[3])