#!/usr/bin/env python3
# coding: utf-8

import sys
import mainFunctions as mf


'''Writes DEGs from DE analysis annotations in text files (.tab)'''


def getGenesAnnotations(filename):

    """
    Reads all genes annotationns from a .txt file and saves it to a dictionary

    Parameter: annotations filename (.tab)

    Return: dictionary with genes as keys and a dictionary with the TAIR (homolog), protein name and GO IDs annotations as value
    """

    genesAnnotations = {}
    
    for line in mf.readFile(filename, True):
        line = line.split('\t')
        genesAnnotations[line[0]] = {'tair': line[1]+'0'.upper(), 'protein': line[3], 'GO': []}
        if line[5] != '':
           genesAnnotations[line[0]]['GO'] = [goid[2:] for goid in line[5].replace('\"', '').split(';')]

    return genesAnnotations


def getDegs4Tables(outdir, allDegs, genesAnnotations):

    """
    Writes gene annotations in tables (.tab) of DEGs from each DE analysis (comparison)

    Parameter:
        directory to output the generated tables
        directory containing all DE results files (.csv)
        annotations filename (.tab)
        
    Return: dictionary with genes as keys and a dictionary with the TAIR (homolog), protein name and GO IDs annotations as value
    """
    
    for comparison in allDegs:

        with open(outdir + comparison + '.txt', 'w') as outfile:
            outfile.write('geneid' +'\t'+ 'tair' +'\t'+ 'protein' +'\t'+ 'GO' +'\n')

            for geneid in allDegs[comparison]:
                outfile.write(geneid +'\t')
                
                # filters genes to write only the annotations from the DEGs of each comparison
                if geneid in genesAnnotations:
                    outfile.write(genesAnnotations[geneid]['tair'] +'\t'+ \
                                  genesAnnotations[geneid]['protein'] +'\t'+ \
                                  ';'.join(genesAnnotations[geneid]['GO']) +'\t')

                else:
                    outfile.write('Unknown gene' +'\t'+ 'Unknown protein' +'\t'+ 'Unknown function' +'\t')

                outfile.write(allDegs[comparison][geneid] + '\n')

        outfile.close()


getDegs4Tables(sys.argv[1], mf.getAllDegsDict(sys.argv[2]), getGenesAnnotations(sys.argv[3]))