#!/usr/bin/env python3
# coding: utf-8

import mainFunctions as mf
import sys


'''Selects the genes with unique homologs, which has higher identity'''


def getAnnotationsIdentities(annotations_filename): 

    """
    Creates a dictionary with transcript ids as keys and a dictionary with gene homologs and identities
    Parameter: tab-delimited annotations file
    Return: dictionary with transcript ids as keys and a dictionary with gene homologs and identities
    """

    trans_annotations = {}

    for line in mf.readFile(annotations_filename, True):
        line = line.strip('\n').split('\t')

        trans_annotations[line[0]] = {'tair': line[1], 'geneName':line[2], 'protName':line[3], 'kegg':line[4], 'goTerms':line[5].replace('\"', ''), 'eVal':line[6], 'identities': line[7]}

    return trans_annotations


def getDegsMaxIsoformTranscriptLength(fasta, allDegsSet):

    """
    Creates a dictionary with all DEGs transcript ID as keys and the respective isoform max length as values
    Parameter: FASTA file with unigenes with the respective max length isoform
    Return: dictionary with all DEGs transcript ID as keys and the respective isoform max length as values
    """

    unigenes_len = {}

    with open(fasta, 'r') as infile:
        sequences = [f.strip('\n') for f in infile.readlines()]

        for i in range(0,len(sequences),2):
            isoform = sequences[i].split(' ')[0][1:]
            transcriptid = '_'.join(isoform.split('_')[:-1])
            lenght = sequences[i].split(' ')[1].split('=')[1]

            if transcriptid in allDegsSet:
                if transcriptid not in unigenes_len:
                    unigenes_len[transcriptid] = lenght

                else:
                    if int(lenght) > int(unigenes_len[transcriptid]):
                        unigenes_len[transcriptid] = lenght

    return unigenes_len


def getTranscriptsHomologs(trans_annotations, unigenes_len):

    """
    Creates a dictionary with TAIRs as values and the most similat (higher identities) transcript IDs as keys
    Parameter: dictionary with transcript ids as keys and a dictionary with gene homologs and identities
    Return: dictionary with TAIRs as values and the most similat (higher identities) transcript IDs as keys
    """

    tair_trans = {}

    for transcriptid in trans_annotations:

        if transcriptid in allDegsSet:

            tair = trans_annotations[transcriptid]['tair']
            identities = trans_annotations[transcriptid]['identities']

            if tair not in tair_trans:
                tair_trans[tair] = {'transcriptid': transcriptid, 'identities': identities, 'len': unigenes_len[transcriptid]}

            if identities == tair_trans[tair]['identities']:
                if int(unigenes_len[transcriptid]) > int(tair_trans[tair]['len']):
                    tair_trans[tair] =  {'transcriptid': transcriptid, 'identities': identities, 'len': unigenes_len[transcriptid]}

            if int(identities) > int(tair_trans[tair]['identities']):
                tair_trans[tair] = {'transcriptid': transcriptid, 'identities': identities, 'len': unigenes_len[transcriptid]}

    trans_tair = {}

    for tair in tair_trans:
        trans_tair[tair_trans[tair]['transcriptid']] = tair

    return trans_tair


def getUniqueTairs(trans_tairs, allDegsDict):

    """
    Creates a dictionary with DEGs with unique annotations from multiple comparisons
    Parameters:
        directory with all DEGs .txt files
        dictionary with transcript ids as keys and a dictionary with gene homologs and identities
        dictionary with TAIRs as values and the most similat (higher identities) transcript IDs as keys
    Return: dictionary with DEGs with unique annotations from multiple comparisons, with transcript ID as keys ands tairs as values
    """

    uniqueDegs = {}

    for comparison in allDegsDict:
        uniqueDegs[comparison] = {}

        for transcriptid in allDegsDict[comparison]:

            if transcriptid in trans_tairs:
                uniqueDegs[comparison][transcriptid] = trans_tair[transcriptid]
    
    return uniqueDegs


def writeUniqueDegs(uniqueDegs, output_dir, trans_annotations, allDegsDict):

    """
    Writes text files with uniquely annotated DEGs annotations
    Parameters:
        dictionary with DEGs with unique annotations from multiple comparisons, with transcript ID as keys ands tairs as values
        directory to output written files
        dictionary with transcript ids as keys and a dictionary with gene homologs and identities
    """

    for comparison in uniqueDegs:

        with open(output_dir + comparison + '.tab', 'w') as outfile:
            outfile.write('Transcript ID' +'\t'+ 'Log2FC' +'\t'+ 'TAIR' +'\t'+ 'Gene Name' +'\t'+ 'Protein Name' +'\t'+ 'KEGG' +'\t'+ 'GO Terms' +'\t'+ 'e-Value' +'\t'+ 'Identities' +'\n')

            for transcriptid in uniqueDegs[comparison]:
                outfile.write(transcriptid +'\t'+ allDegsDict[comparison][transcriptid] +'\t')

                for attribute in trans_annotations[transcriptid]:
                    outfile.write(trans_annotations[transcriptid][attribute] +'\t')
                outfile.write('\n')


# DEGs directory
allDegsDict = mf.getAllDegsDict(sys.argv[1], '\t', -1)
allDegsSet = mf.getAllDegsSet(sys.argv[1], '\t')

# annotations file
trans_annotations = getAnnotationsIdentities(sys.argv[2])

# transcriptome file (FASTA)
trans_tair = getTranscriptsHomologs(trans_annotations, getDegsMaxIsoformTranscriptLength(sys.argv[3], allDegsSet))

uniqueDegs = getUniqueTairs(trans_tair, allDegsDict)

# output directory
writeUniqueDegs(uniqueDegs, sys.argv[4], trans_annotations, allDegsDict)