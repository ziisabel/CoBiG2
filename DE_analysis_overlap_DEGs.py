#!/usr/bin/env python3
# coding: utf-8

import sys
import mainFunctions as mf


'''Writes a tab-delimited text file (.tab) with the overlapping results from DESeq2 and edgeR'''


def getAllDegsFromDE(filename, sep):

    """
    Creates a dictionary with the DE results from DESeq2 or edgeR

    Parameter: DE results file from DESeq2 or edgeR
    
    Return: dictionary with DEGs IDs as keys and log2FC and FDR as values
    """

    allDegs = {}
    for line in mf.readFile(filename, True):
        line = line.replace('\"', '').strip('\n').split(sep)
        allDegs[line[0]] = {'log2FC': line[3], 'FDR': line[-1]}

    return allDegs


def getOverlappingDEGs(allDegsDESeq2, allDegsEdgeR):

    """
    Creates a dictionary with DEGs detected simultaneously by DESeq2 and edgeR

    Parameter:
        output from the function getAllDegsFromDESeq2()
        output from the function getAllDegsFromEdgeR()

    Return: dictionary with DEGs IDs as keys and log2FC as values
    """

    overlappingDegs = {}

    for geneid in allDegsDESeq2:
        if geneid in allDegsEdgeR:
            overlappingDegs[geneid] = allDegsDESeq2[geneid]['log2FC']

    return overlappingDegs


def writeOverlappingDegs(overlappingDegs, filename):

    """
    Writes a text file (.tab) with the overlapping DE results from DESe2 and edgeR

    Parameter:
         output from the function getOverlappingDEGs()
        output filename with overlapping DE results from DESe2 and edgeR
    
    Output: .tab file with the overlapping DE results from DESe2 and edgeR
    """
    
    with open(filename, 'w') as outfile:
        outfile.write('geneID\tlog2FC')

        for geneid in overlappingDegs:
            outfile.write(geneid +'\t'+ overlappingDegs[geneid])

    outfile.close()


allDegsDESeq2 = getAllDegsFromDE(sys.argv[1])
allDegsEdgeR = getAllDegsFromDE(sys.argv[2])

writeOverlappingDegs(getOverlappingDEGs(allDegsDESeq2, allDegsEdgeR), sys.argv[3])