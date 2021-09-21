#!/usr/bin/env python3
# coding: utf-8

import csv
import os
import os.path
import sys
import mainFunctions as mf


'''Formats the .cvs files retrieved from gProfiler ORA for easier readbility and extraction for REVIGO and/plotting'''


def getGProfiler(file):
    """
    Creates a dictionary with the enrichment results retrieved from gProfiler ORA

    Parameter: .csv file retrieved from gProfiler ORA

    Return: dictionary with the enrichment results for each GO term of the input file (keys)
    """

    gProfiler = {}

    # saves file content into a dictionary with GO terms as keys and the enrichment information as value
    for line in mf.readFile(file, True):
        line = [l.replace('\"', '') for l in line.split(',')]
        for l in line[1:]:
            if 'GO:' in l or 'KEGG:' in l:
                index = line.index(l)
                gProfiler[line[index]] = {'Database': line[0],
                                        'Description': ('').join(line[1:index]),
                                        'Term ID': line[index], 'padj': line[index + 1],
                                        'Counts': line[index + 5],
                                        'GeneRatio': str(float(line[index + 5]) / float(line[index + 3])),
                                        'Comparison': file.strip('.csv')}
    
    return gProfiler


def cleanGProfiler(path, outdir_cleaned, outdir_revigo):
    """
    Writes gProfiler dictionaries' contents and REVIGO input into .tab files

    Parameters: 
        path to gProfiler results (.csv) directory
        path to gProfiler output directory
        path to REVIGO input files

    Output:
        gProfiler output (cleaned) results in .tab format
        REVIGO input files in .tab format
    """

    gProfiler = []

    for file in [f for f in os.listdir(path) if not f.startswith('.')]:

        gProfiler = getGProfiler(path + file)

        # writes gProfiler dictionaries' contents into .tab files
        with open(outdir_cleaned + file.strip('.csv') + '.tab', 'w') as outfile:
            writer = csv.writer(outfile, dialect="excel-tab")
            writer.writerow(['Database', 'Description', 'Term ID', 'padj', 'Counts', 'GeneRatio', 'Comparison'])

            for termID in gProfiler:
                writer.writerow(gProfiler[termID].values())
        outfile.close()

        # writes REVIGO input into .tab files
        with open(outdir_revigo + file.strip('.csv') + '.tab', 'w') as outfile:
            writer = csv.writer(outfile, dialect="excel-tab")
            writer.writerow(['Term ID', 'padj'])

            for termID in gProfiler:
                writer.writerow([gProfiler[termID]['Term ID'], gProfiler[termID]['padj']])
        outfile.close()


cleanGProfiler(sys.argv[1], sys.argv[2], sys.argv[3])