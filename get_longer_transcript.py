#!/usr/bin/env python3
# coding: utf-8

import sys
import mainFunctions as mf


'''Creates a sequences file (.fasta) containing only the longer transcripts sequences of each gene in a transcriptome'''


def getTairSequences(filename, degsTairs):

    """
    Creates a dictionary with DEGs' homologs TAIRs as keys and their longer sequences according to the transcriptome as values

    Parameters:
        Transcriptome .fasta file
        dictionary with genes IDs as keys and the respective TAIRs as values

    Return: dictionary with DEGs' homologs TAIRs as keys and their longer sequences according to the transcriptome as values
    """

    sequencesLen = {}
    transcriptsSequences = {}
    
    # read transcriptome .fasta file and saves it into a list
    rawSequences = mf.readFile(filename, False)
    for line in rawSequences:
        if '>' in line:
            lineList = line.split(' ')
            # find the gene IDs in the transcriptome
            geneid = lineList[0][1:-3]
            # check which genes are DEGs
            if geneid in degsTairs:
                # find the transcript IDs for each gene in the transcriptome
                transcriptid = lineList[0][1:]
                if geneid not in sequencesLen:
                    sequencesLen[geneid] = {}
                # save all the transcripts sequence lenghts of a gene into a dictionary
                sequencesLen[geneid][transcriptid] = int(lineList[1][4:])
                # save all transcripts sequence into a dictionary
                transcriptsSequences[transcriptid] = rawSequences[rawSequences.index(line)+1]

    degsTairsSequences = {}
    for geneid in sequencesLen:
        # for each DEG, find which respective transcript sequence is longer and saves it into a dictionary with homologs TAIRs as keys
        degsTairsSequences[degsTairs[max(sequencesLen[geneid], key=lambda key: sequencesLen[geneid][key])[:-3]]] = transcriptsSequences[transcriptid]
    
    return degsTairsSequences


def writeTairsFasta(filename, degsTairsSequences):

    """
    Writes a .fasta file containing only the longer transcripts sequences of each gene in the transcriptome

    Parameters:
        .fasta filename to output
        dictionary with DEGs' homologs TAIRs as keys and their longer sequences according to the transcriptome as values

    Output: .fasta file containing only the longer transcripts sequences of each gene in the transcriptome
    """

    with open(filename, 'w') as outfile:
        for tair in degsTairsSequences:
            outfile.write('>' + tair +'\n'+ degsTairsSequences[tair] +'\n')

        outfile.close()


degsTairsSequences = getTairSequences(sys.argv[1], mf.getDegsTairs(sys.argv[2]))
writeTairsFasta(sys.argv[3], degsTairsSequences)