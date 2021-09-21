#!/usr/bin/env python3
# coding: utf-8

import sys


'''Converts a multiple-line per gene .fasta file into a single-line .fasta file'''


def writeCleanFastaGenome(file_in, file_out, pattern):
    """
    Converts a multiple-line per gene fasta file to single-line

    Parameters:
        multiple-line per gene fasta file to convert
        single-line per gene fasta file to write
        character or string before gene ID

    Output: single-line per gene fasta file
    """

    # reads multiple-line per gene fasta file and saves to a list of file line
    with open(file_in) as infile:
    	rawSequences = infile.readlines()
    infile.close()

    with open(file_out, 'w') as outfile:

        for line in rawSequences:
            # finds the lines with gene ID
            if pattern in line:
                linelist = line.split(' ')
                
                # writes gene ID on the output file in the following format: > geneID
                for attribute in linelist:
                    if pattern in attribute:
                        geneid = '>' + attribute.replace(pattern, '')
                        outfile.write(geneid + '\n')
                        break
            
                # finds and joins multiple-line sequences per gene into a single line and writes it to the output file bellow the respective gene ID
                for i in range(rawSequences.index(line) + 1, len(rawSequences)):
                    if rawSequences[i].startswith('>'):
                        fragments = [fragment.strip('\n') for fragment in rawSequences[rawSequences.index(line)+1:i-1]]
                        sequence = ('').join(fragments)
                        outfile.write(sequence + '\n')
                        break
                    elif i == len(rawSequences)-1:
                        for j in range(i, 0, -1):
                            if rawSequences[j].startswith('>'):
                                fragments = [fragment.strip('\n') for fragment in rawSequences[j+1:len(rawSequences)-1]]
                                sequence = ('').join(fragments)
                                outfile.write(sequence + '\n')
                                break

        outfile.close()


writeCleanFastaGenome(sys.argv[1], sys.argv[2], sys.argv[3])