#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import mainFunctions as mf


'''Writes text annotation files (.tab or .txt) with gene IDs, ORF names, gene names, protein names and GO terms'''


def getGeneName(tid):

	""""
	Searches Uniprot for the gene name of an OFR name

	Parameter: ORF name

	Return: Gene name
	"""
	
	geneName = tid
	
	try:

		# search uniprot for the ORF name
		htmlResponse_info = mf.urlOpenWithRetry('https://www.uniprot.org/uniprot/?query=' + tid + ',&columns=id,entry%20name,go&format=tab')
		stringLines_info = htmlResponse_info.data.decode('utf-8').split('\n')
		
		# get uniprot acession number of the ORF name
		if stringLines_info != []:
			for line in stringLines_info[1:]:
				acession = line.split('\t')[0]
				break
		
			# search uniprot for the respective uniprot acession number annotation
			htmlResponse_info = mf.urlOpenWithRetry('https://www.uniprot.org/uniprot/' + acession + '.txt')
			stringLines_info = htmlResponse_info.data.decode('utf-8').split('\n')

			if stringLines_info != []:
				for line in stringLines_info:
				
					# get the gene name associated with the uniprot acession number
					if line.startswith('GN'):
						geneName = line.split('{')[0].split('=')[1].strip(' ')
						break

		return geneName
	
	except:
		return geneName
			


def getGenesProtNames(filename):

	""""
	Searches an annotations file (.gff3) for the protein names of each ORF names and saves them into a dictionary
	
	Parameter: annotations filename (.gff3)
	
	Return: dictionary with ORF names as keys and protein names as values
	"""

	protNames = {}

	# search protein names in a list of .gff3 annotations and saves them as values of a dictionary with ORF names as keys
	for line in mf.readFile(filename, True):
		if 'polypeptide' in line:
			line = line.strip('\n').split('\t')[-1]
			protName = line[line.index('PRODUCT:'):].split(',')[0].split(':')[1].replace('%3B', ';').replace('%2C', ',')
			tid = line[line.index('note='):].split('~')[0].split('=')[1]
			
			protNames[tid] = protName
	
	return protNames


def getGenesAnnotations(genesIds_file, protNames):

	""""
	Creates a dictionary with gene ID as keys and ORF name, gene name, protein name and GO terms as values
	
	Parameter:
		.txt file with gene ID and the respective ORF name (one pair per line)
		dictionary with ORF names as keys and the respective protein names as values
	
	Return: dictionary with gene ID as keys and ORF name, gene name, protein name and GO terms as values
	"""

	genesInfo = {}

	# reads a .txt file with gene ID and the respective ORF name (one pair per line)
	for line in mf.readFile(genesIds_file, False):
		
		# get each gene ID and respective ORF name
		line = line.split('\t')
		gene = line[1]
		tid = line[0]
		
		# create a dictionary with gene ID as keys and ORF name, gene name, protein name and GO terms as values
		genesInfo[gene] = {'tid': '', 'geneName': '', 'protName': '', 'goTerms': ''}
		genesInfo[gene]['tid'] = tid
		genesInfo[gene]['geneName'] = getGeneName(tid.replace('T', '_T'))
		genesInfo[gene]['protName'] = protNames[tid]

		# search the annotations file (.gff3) for GO terms and adds them to the dictionary
		goTerms = os.popen('cat STAR/annotations/c_canephora_gene_structural_and_functional_annotation.gff3 | grep ' + gene + ' | grep -o GO:[0-9]*').readlines()
		goTerms = [goTerm.strip('\n') for goTerm in goTerms]
		goTerms = (';').join(goTerms)
		genesInfo[gene]['goTerms'] = goTerms
	
	return genesInfo
	

def writeGenesAnnotations(annotations_file, genesInfo):

	""""
	Writes a text file (.tab or .txt) with gene IDs, ORF names, gene names, protein names and GO terms
	
	Parameter:
		.txt file with tab separated genes annotations (one per line)
		dictionary with gene ID as keys and ORF name, gene name, protein name and GO terms as values
	
	Output: .tab file with gene IDs, ORF names, gene names, protein names and GO terms (1 gene per line)
	"""

	with open(annotations_file, 'w') as outfile:
		# write output file header (1st line)
		outfile.write('ID\tTID\tGene_Name\tProtein_Name\tGO_Terms\n')

		for gene in genesInfo:

			# create a list with the gene values from the dictionary
			row = [gene]
			for attribute in genesInfo[gene]:
				row.append(genesInfo[gene][attribute])

			# writes the list as a tab separated line in the output file (.tab)
			row = ('\t').join(row)
			outfile.write(row + '\n')

	outfile.close()


writeGenesAnnotations(sys.argv[1], getGenesAnnotations(sys.argv[2], getGenesProtNames(sys.argv[3])))