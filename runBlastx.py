#!/usr/bin/env python
# coding: utf-8

import os
import sys
import mainFunctions as mf


'''Runs blastx for a set of genes in a reference genome, mapping them to a specific database (e.g. swissprot)'''
    

def runBlastxFromFasta(fastas_outdir, genome_fasta, allDegs, blastx_outdir, database_dir):

    """
    Parses a genome into multiple fasta files (one per gene) and blastx them into a specified database

    Parameters:
        directory to temporarly output fasta files
        reference genome (.fasta)
        set of DEGs output from the function getAllDegs()
        directory to output the blastx files
        directory containing the specified database

    Output: multiple blastx files (one per gene)
    """

    genome = {}

    # save all genome sequences into a list
    with open(genome_fasta, 'r') as infile:
        sequences = [f.strip('\n') for f in infile.readlines()]

        for i in range(0,len(sequences),2):
            # get all genes of the genome
            gene = sequences[i].strip('>')
            # save the longest transcript of each gene as a value of the genome dictionary (genes as keys)
            if gene not in genome:
                genome[gene] = '>' + gene + '\n' + sequences[i+1]
            else:
                if len(sequences[i+1]) > len(genome[gene]):
                    genome[gene] = '>' + gene + '\n' + sequences[i+1]

    for gene in genome:
        if gene in allDegs:

            # write a .fasta file for each gene containing its respective sequence (longest transcript)
            fastaFilename = gene + '.fasta'
            with open(fastas_outdir + fastaFilename, 'w') as outfile:
                outfile.write(genome[gene])
            outfile.close()

            blastxFilename = gene + '.blastx'

            # show progress
            print(fastaFilename, '---->', blastxFilename)

            # run blastx for each gene fasta file and save the result in a .blastx file
            os.system('blastx -query ' + fastas_outdir + fastaFilename + ' -out ' + blastx_outdir + blastxFilename + ' -db ' + database_dir + ' -evalue 1000')

            # remove the already used .fasta file to save space
            os.system('rm ' + fastas_outdir + fastaFilename)


runBlastxFromFasta(sys.argv[1], sys.argv[2], mf.getAllDegsSet(sys.argv[3]), sys.argv[4], sys.argv[5])

print('Total blastx files:', len(os.listdir(sys.argv[4])))
print('DONE')