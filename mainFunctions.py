#!/usr/bin/env python3
# coding: utf-8

import os
from urllib3 import Retry, PoolManager


'''Set of functions that are called in multiple scripts'''


def readFile(filename, header):
    """
    Reads input file and saves content into a list

    Parameter:
        full file name with path (.txt, .tab or .csv)
        boolean to indicate if the first line is the header (True or False)

    Return: list with one file line per item
    """

    with open(filename, 'r') as infile:
        if header:
            infile.readline()
        content = [i.strip('\n') for i in infile.readlines()]
    infile.close()

    return content


def getAllDegsDict(indir, sep, logFC_index):

    """
    Reads DEGs from all DE results files and saves them as a dictionary with comparisons (DE) as keys

    Parameter:
        directory containing all DE results files (.csv)
        file separator (single character; e.g: comma, slash, etc.)
        index of the log2FC value in each line of the file (integer)

    Return: set with all DEGs from all comparison
    """
    
    allDegs = {}

    for filename in [f for f in os.listdir(indir) if not f.startswith('.')]:

        degs = {}

        for line in readFile(indir + filename, True):
            line = line.strip('\n').split(sep)
            # add an entry for each gene to the degs dictionary {gene: log2FC}
            degs[line[0]] = line[logFC_index]
        
        # add the degs dictionary as value to the comparison dictionary {comparison: {degs}}
        allDegs[filename.strip('.csv')] = degs

    return allDegs


def getAllDegsSet(indir):
    
    """
    Reads all DEGs from DE results files and saves them as a unique set

    Parameter: directory containing all DE results files (.csv)

    Return: set with all DEGs from all comparison
    """
    
    allDegs = set()

    for filename in [f for f in os.listdir(indir) if not f.startswith('.')]:
        for line in readFile(indir + filename, True):
            allDegs.add(line.strip('\n').split('\t')[0])

    return allDegs


def getDegsTair(degsTair_file):

    """
    Generates a dictionary with DEGs as keys and the respective TAIRs homologs (mapped through blast) as values

    Parameter: tab-delimited file (.txt or .tab) with DEGs and the respective TAIRs homologs (1 pair per line)

    Return: dictionary with DEGs as keys and the respective TAIRs homologs as values
    """

    degsTairs = {}

    for line in readFile(degsTair_file, True):
        line = line.split('\t')
        degsTairs[line[0]] = line[1]

    return degsTairs


def urlOpenWithRetry(url):

    """
    Make a request to access a resource on the server

    Parameter: Uniform Resource Locator (URL)

    Return: Response object from the server (http.request)
    """

    http = PoolManager()
    response = http.request('GET', url, retries=Retry(3))

    return response