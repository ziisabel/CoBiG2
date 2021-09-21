#!/usr/bin/env python3
# coding: utf-8

import sys
import json
import mainFunctions as mf


'''Generates all the child terms of specific GO IDs from a set of annotated genes'''


def getTermChildren(goid, goChildren):

    """
    Fills an empty set with all the child terms of a specific GO ID, recursively

    Parameter:
        GO ID
        empty set to be filled with child terms, recursively
    """

    try:

        goChildren.add(goid)

        # Search the child terms of a specific GO ID using Quick GO API
        htmlResponse_info = mf.urlOpenWithRetry('https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/'+ goid.replace(':', '%3A') + '/children')
        stringLines_info = htmlResponse_info.data.decode('utf-8').split('\n')

        results = json.loads(stringLines_info[0])['results'][0]
        if 'children' in results:
            children = results['children']
            
            # recursively add the child terms of each child term found previously
            for child in children:
                goChildren.add(child['id'])
                getTermChildren(child['id'], goChildren)

        return 

    except:
        return


def getGenesAnnotations(filename, sep):

    """
    Reads all genes annotationns from a .txt file and saves it to a dictionary

    Parameter:
        annotations filename (.txt)
        file separator (single character; e.g: comma, slash, etc.)

    Return: dictionary with genes as keys and a dictionary with the gene name, protein name and GO IDs annotations as value
    """

    genesAnnotations = {}
    
    for line in mf.readFile(filename, True):
        line = line.split(sep)
        genesAnnotations[line[0]] = {'geneName':line[2], 'protein': line[3]}
        if line[4] == '':
            genesAnnotations[line[0]]['GO'] = []
        else:
            genesAnnotations[line[0]]['GO'] = line[4].split(';')

    return genesAnnotations


def getspecificGOdegs(allDegs, genesAnnotations, goChildren):

    """
    Creates a dictionary with gene IDs as keys and all the child terms of their respective GO IDs as values

    Parameter:
        dictionary with genes as keys and a dictionary with the gene name, protein name and GO IDs annotations as value
        set with all the child terms of a specific GO ID

    Return: dictionary with gene IDs as keys and all their GO IDs (direct and child) as values
    """

    specificGOdegs = {}

    for comparison in allDegs:
        for geneid in allDegs[comparison]:
            for goid in goChildren:
                if goid in genesAnnotations[geneid]['GO']:
                    specificGOdegs[geneid] = genesAnnotations[geneid]
                    specificGOdegs[geneid]['GO'] = goid

    return specificGOdegs


def writeSpecificGODegs(filename, allDegs, specificGOdegs, genes_tair):

    """
    Writes a tab-delimited text file (.txt) with comparisons (DE), gene IDs, protein names, GO IDs and, if any, TAIRs homologs
    
    Parameter:
        output filename (.txt)
        output of the function getAllDegs()
        output of the function getspecificGOdegs()
        output of the function getGenesTairs()

    Output:
    """

    with open(filename, 'w') as outfile:
        
        for comparison in allDegs:
            for geneid in allDegs[comparison]:
                if geneid in specificGOdegs:
                    outfile.write(comparison +'\t'+ geneid +'\t')

                    # if available, write the gene TAIR homolog
                    if geneid in genes_tair:
                        outfile.write(genes_tair[geneid] +'\t')
                    else:
                        outfile.write('' +'\t')
                    outfile.write(specificGOdegs[geneid]['protein'] +'\t'+ specificGOdegs[geneid]['GO']+'\t'+ allDegs[comparison][geneid] +'\n')

    outfile.close()

allDegs = mf.getAllDegsDict(sys.argv[1], sys.argv[2], sys.argv[3])

children = set()
for goid in sys.argv[4]:
    getTermChildren(goid, children)

goTerms_degs = getspecificGOdegs(allDegs, getGenesAnnotations(sys.argv[5]), children)
writeSpecificGODegs(sys.argv[6], allDegs, goTerms_degs, mf.getDegsTair(sys.argv[7]))