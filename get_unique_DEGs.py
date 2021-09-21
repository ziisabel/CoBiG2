#!/usr/bin/env python3
# coding: utf-8


'''Selects the longer transcript of each given gene'''


def getUniqueDegs(degGenes):

    """
    Finds transcripts mapped to the same gene (through blast) and filters them by sequence lenght, selecting the longer

    Parameters: dictionary with transcript ids as keys and protein names and indentities percentages from blast as values.

    Return: dictionary of transcripts uniquely annotated to different genes
    """

    # creates a dictionary with gene annotations (protein names) as keys and the list of respective transcripts as values
    protNamesDegs = {}
    for gene in degGenes:
        
        if degGenes[gene]['protName'] not in protNamesDegs:
            protNamesDegs[degGenes[gene]['protName']] = [gene]
        else:
            protNamesDegs[degGenes[gene]['protName']].append(gene)
    
    # creates a dictionary with gene annotations (protein names) as keys and a list of the transcripts' identities percentages from blast as values
    protIdentities = {}
    for protName in protNamesDegs:
        identitiesList = []
        for gene in protNamesDegs[protName]:
            identitiesList.append(int(degGenes[gene]['identities']))
        protIdentities[protName] = identitiesList

    # creates a dictionary with gene annotations (protein names) as keys and the respective transcript with higher identity as values
    uniqueDegsDict = {}
    for protName in protIdentities:
        if len(protIdentities[protName]) > 1:
            maxIndex = protIdentities[protName].index(max(int(sub) for sub in protIdentities[protName]))
            uniqueDegsDict[protName] = protNamesDegs[protName][maxIndex]
        else:
            uniqueDegsDict[protName] = protNamesDegs[protName][0]

    # creates a dictionary with filtered transcripts as keys and the respective annotations (protein names) as values
    uniqueDegs = {}
    for protName in uniqueDegsDict:
        uniqueDegs[uniqueDegsDict[protName]] = degGenes[uniqueDegsDict[protName]]
        
    return uniqueDegs