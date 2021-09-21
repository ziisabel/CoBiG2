#!/usr/bin/env python3
# coding: utf-8

import mainFunctions as mf
import sys

'''Generates a text file (.tab) with a DE comparison name, transcript/gene ID, Arabidopsis thaliana homolog TAIR,
the respective TF family name retrieved from iTAK (Plant Transcription factor & Protein Kinase Identifier and Classifier) and the type of regulation'''


def getTairsTF(filename):

    """
    Creates a dictionay with TAIR as keys as the respective TF family as values

    Parameter: filename containing a list of TAIRs with all the respective TFs family names (one gene per line)

    Return: dictionay with TAIR as keys as the respective TF family as values
    """

    with open(filename, 'r') as infile:
        content = [line.strip('\n') for line in infile.readlines()]
    infile.close()

    tairsTF = {}

    for line in content:
        if '.' in line:
            line = line.split('\t')
            tairsTF[line[0][:-2]] = line[1]
        else:
            print(line)

    return tairsTF


def getTairDegs(filename, allDegsSet):

    """
    Creates a dictionary with a TAIRs as keys and the respective gene/transcript IDs as values

    Paratemer: .txt annotations file containing at least gene/transcript IDs and 

    Return: dictionary with a TAIRs as keys and the respective gene/transcript IDs as values
    """

    tairDegs = {}
    
    for line in mf.readFiles(filename, True):
        line = line.split('\t')
        tair = line[1].upper()
        if tair in allDegsSet:
            tairDegs[tair] = line[0]
        
    return tairDegs


def getDegsTF(tairsTF, tairDegs):

    """
    Creates a dictionary with DEGs ID as keys and another dictionary with the respective TAIR and TF as values

    Parameter:
        output of the getTairsTF() function
        output of the getTairDegs() function
    
    Return: dictionary with DEGs ID as keys and another dictionary with the respective TAIR and TF as values
    """

    degsTF = {}

    for tair in tairsTF:
        if tair in tairDegs:
            degsTF[tairDegs[tair]] = {'tair': tair, 'TF': tairsTF[tair]}
    
    return degsTF


def writeDegsTF(filename, degsTF):

    """
    Creates a .txt tab-delimited file with gene IDs, TAIRs and the respective TFs (one gene per line)
    Parameter:
        fullname of the file to output (.txt)
        output of the getDegsTF() function
    """

    with open(filename, 'w') as outfile:
        for geneid in degsTF:
            outfile.write(geneid +'\t'+ degsTF[geneid]['tair'] +'\t'+ degsTF[geneid]['TF'] +'\n')


def getTairsTFregulation(allDegs, degsTF):

    """
    Creates a dictionary with multiple DE comparisons, with gene IDs as keys and a dicionary with TF families and regulation type as values
    
    Parameter:
        dictionary with comparisons (DE) as keys and a dictionary as value with gene ID as keys and log2FC as values
        output of the getDegsTF() function
    
    Return: dictionary with multiple DE comparisons, with gene IDs as keys and a dicionary with TF families and regulation type as values
    """

    tairsTFregulation = {}

    for comparison in allDegs:
        tairsTFregulation[comparison] = {}
        for geneid in allDegs[comparison]:
            if geneid in degsTF:
                tairsTFregulation[comparison][geneid] = degsTF[geneid]
                if float(allDegs[comparison][geneid]) > 0:
                    tairsTFregulation[comparison][geneid]['regulation'] = 'up'
                else:
                    tairsTFregulation[comparison][geneid]['regulation'] = 'down'

    return tairsTFregulation


def writeTairsTFregulation(filename, tairsTFregulation):

    """
    Creates a .txt tab-delimited file with the getTairsTFregulation() content (one gene per line)
   
    Parameter:
        fullname of the file to output (.txt)
        output of the getTairsTFregulation() function
    """

    with open(filename, 'w') as outfile:
        for comparison in tairsTFregulation:
            for geneid in tairsTFregulation[comparison]:
                row = [comparison, geneid]
                for attribute in tairsTFregulation[comparison][geneid]:
                    row.append(tairsTFregulation[comparison][geneid][attribute])
                outfile.write('\t'.join(row) + '\n')


tairs_TF = getTairsTF(sys.argv[1])
annotations = getTairDegs(sys.argv[2], mf.getAllDegsSet(sys.argv[3], '\t'))
degsTF = getDegsTF(tairs_TF, annotations)
writeDegsTF(sys.argv[3], degsTF)

allDegs = mf.getAllDegsDict(sys.argv[3], '\t', 3, '.tab')
tairsTFregulation = getTairsTFregulation(allDegs, degsTF)
writeTairsTFregulation(sys.argv[4], tairsTFregulation)