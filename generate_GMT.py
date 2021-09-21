#!/usr/bin/env python3
# coding: utf-8

import csv
import os
import sys
import mainFunctions as mf


'''Generates personalized GMT files to use as input in gProfiler to perform over-expression enrichment analysis (ORA)'''


def getExpressedGenes(counts_dir):

    """
    Creates a dictionary with gene expression by sample

    Parameter: directory with all count files (.csv, .tab or .txt)

    Return: a dictionary with samples as keys and a list of the respective expressed genes as values
    """

    # creates a dictionary with samples as key and an empty set as value
    allExpressedGenes = {}
    files = [file for file in os.listdir(counts_dir)]

    for file in files:
        sample = file.split('.')[0]
        if sample not in allExpressedGenes:
            allExpressedGenes[sample] = set()

        # checks if each gene count is above zero and if so adds it to the respective sample set on the dictionary
        for line in mf.readFile(counts_dir + file, True):
            line = line.split('\t')
            if int(line[1]) > 0:
                allExpressedGenes[sample].add(line[0])

    return allExpressedGenes


def createGMT(filename_in, samplesExpressedGenes, filename_out):

    """
    Generates a personalized GMT file for over-expression enrichment analysis (ORA)

    Parameter:
        GMT file with all genes
        list of genes expressed by treatment and control
        filtered GMT filename

    Output: .tab file with the filtered GMT (treatment + control)
    """

    # creates a dictionary with GO IDs as keys and GO description and genes annotated by it as values
    allGmt = {}
    for line in mf.readFile(filename_in, False):
        line = line.split('\t')
        allGmt[line[0]] = {'goName': line[1], 'genes': line[2:]}

    # creates a dictionary representing a subset of allGmt with only genes expressed by treatment and control samples
    expressedGMT = {}
    for goId in allGmt:
        expressedGMT[goId] = {'goName': allGmt[goId]['goName'], 'genes': []}

        for gene in allGmt[goId]['genes']:
            if gene in samplesExpressedGenes:
                expressedGMT[goId]['genes'].append(gene)

    # writes a .tab file with the filtered GMT information (treatment + control)
    with open(filename_out, 'w') as outfile:
        writer = csv.writer(outfile, dialect="excel-tab")

        for goId in expressedGMT:
            if expressedGMT[goId]['genes'] != {}:
                row = [goId, expressedGMT[goId]['goName']]

                for gene in expressedGMT[goId]['genes']:
                    row.append(gene)
                writer.writerow(row)
    
    outfile.close()


def generateAllGmtFiles(allExpressedGenes, treatment, control, gmt_file, output_path):

    """
    Runs the createGMT() function for a pair of specific treatment and control samples

    Parameters:
        dictionary of with gene expression by sample
        treatment sample
        control sample
        GMT file with all genes
        path to write output
    """

    # creates a dictionary with the filtered gene expression of the treatment and control samples
    samplesExpressedGenes = allExpressedGenes[treatment]
    samplesExpressedGenes.update(allExpressedGenes[control])

    # call the createGMT() funtion for the selected treatment and control samples
    createGMT(gmt_file, samplesExpressedGenes, output_path + treatment + 'vs' + control + '.gmt')


generateAllGmtFiles(getExpressedGenes(sys.argv[1]), sys.argv[2], sys.argv[3])