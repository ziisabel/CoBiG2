#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import mainFunctions as mf


'''Creates text files (.tab) to use as input in gProfiler to perform enrichment analysis (ORA)'''


def getAllDegsIds(DE_directory):

    """
    Creates a dictionary for all the comparisons (DE) with the respective DEGs and their log2FC values, grouped by type of regulation
    
    Parameter: directory with all comparison (DE) files
    
    Return: dictionary with all the comparisons (DE), the respective DEGs and their log2FC values, grouped by type of regulation
    """

    allDegsIds = {}
    regulation = ['up', 'down']
    # create a list with all comparison (DE) filenames
    files = [file for file in os.listdir(DE_directory) if file[0].isdigit()]

    for file in files:
        comparison = file.split('_')[0]
        # create a dictionary with comparisons (DE) as values and regulation type (up or down) as values
        allDegsIds[comparison] = {reg: {} for reg in regulation}

        # read each comparioson (DE) file
        for line in mf.readFile(DE_directory + file, True):
            line = line.split('\t')

            # check the log2FC of each gene in a comparison and save them in the respective regulation key of the dictionary
            if float(line[1]) > 0:
                allDegsIds[comparison]['up'][line[0]] = line[1]
            else:
                allDegsIds[comparison]['down'][line[0]] = line[1]
        
        for reg in regulation:
            # orderes the regulation dictionaries from the most to the least DE genes (descending positive values and ascending negative values)
            allDegsIds[comparison][reg] = {k: v for k, v in sorted(allDegsIds[comparison][reg].items(), key=lambda item: item[1], reverse=True)}

    return allDegsIds


def writeOrderedDegs(allDegsIds, gProfiler_directory):

    """
    Creates the input files .tab to perform enrichment analysis (ORA) in gProfiler
    
    Parameter:
        dictionary with all the comparisons (DE), the respective DEGs and their log2FC values, grouped by type of regulation
        directory to save the output files
    
    Output: input files to perform enrichment analysis (ORA) in gProfiler
    """

    for comparison in allDegsIds:
        for regulation in allDegsIds[comparison]:
            # assign the output filename for each comparison (DE)
            filename = gProfiler_directory + comparison + '_' + regulation + '.tab'
            
            # create a .tab file for each comparison (DE)
            with open(filename, 'w') as outfile:
                
                # write the input information for gProfiler in each file
                for deg in allDegsIds[comparison][regulation]:
                    outfile.write(deg +'\t'+ allDegsIds[comparison][regulation][deg])

            outfile.close()


writeOrderedDegs(getAllDegsIds(sys.argv[1]), sys.argv[2])