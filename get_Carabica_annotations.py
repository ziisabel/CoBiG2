#!/usr/bin/env python3
# coding: utf-8

import sys
import mainFunctions as mf


'''Generates a text file (.txt) with gene name and the respective protein annotation from the C. arabica reference genome'''


def getAnnotations(filename):

    """
    Creates a dictionary with genes as keys and proteins as values extracted from the C. arabica reference genome
    Paramenter: GCF_003713225.1_Cara_1.0_cds_from_genomic.fna (C. arabica genome annotations)
    Return: dictionary with genes as keys and proteins as values from C. arabica
    """

    annotations = {}

    # search for genes and protein names and saves them into a dictionary with genes as keys and proteins as values
    for line in mf.readFile(filename, False):
        if line.startswith('>'):
            geneid = line.split('gene=')[1].split(']')[0]
            protein = line.split('protein=')[1].split(']')[0]
            annotations[geneid] = protein

    return annotations


def writeAnnotations(filename, annotations):

    """
    Creates a text file (.txt) with gene and protein names from the annotations dictionary
    Parameter: output filename (.txt)
    Output: .txt file with gene and protein names from the C. arabica reference genome
    """

    with open(filename, 'w') as outfile:
        # for each key (gene) in the dictionary writes a line with the gene and protein annotation, tab-delimited
        for geneid in annotations:
            outfile.write(geneid + '\t' + annotations[geneid] + '\n')
    outfile.close()


annotations = getAnnotations(sys.argv[1])
writeAnnotations(sys.argv[2], annotations)