#!/usr/bin/env python3
# coding: utf-8

import sys
import os 
import os.path
import mainFunctions as mf


'''Writes a tab-delimited text file (.tab) with Uniprot annotations for blastx files of a set of DEGs'''


def getBlastxValues(content, i):

    """
    Creates a dictionary with the e-value and identity values for a given blastx result

    Parameters:
        list containing the description lines of a BLASTX file (1 line per item)
        index of the content list that contains an organism name (integer)

    Return: dictionary with e-value and identities values for a given blastx result
    """

    for j in range(i, len(content)):
        if content[j].startswith(' Score'):
            evalue = content[j].split(',')[1].split(' ')[-1]
            identities = content[j+1].split(',')[0].split(' ')[-1][1:-2]
            break
            
    return {'evalue': evalue, 'identities': identities}


def getUniprotHits(filename, species):

    """
    Reads a .blastx file and searches for the first entry of a specified species that passes the e-value (<0.001) and identity (>40%) thresholds

    Parameter:
        .blastx filename to read
        name of the species to search in the .blastx file

    Return: a dictionary with an uniprot acession number as key and the respective e-value and identity as values
    """

    uniprot_hits = {}
    if not os.path.isfile(filename):
        print(filename, "doesn't exist!")
        return uniprot_hits # None
    
    # get the description lines of a BLASTX file
    with open(filename, 'r') as infile:
        header = infile.readline()
        content = infile.readlines()
    infile.close()

    # sanity check that this is a BLASTX output text file
    if not header.startswith("BLASTX"):
        #raise TypeError("Not a BLASTX file")
        print("Not a BLASTX file")
        return None
    
    organism = ''

    for line in content:

        if line.startswith(">"):
            for i in range(content.index(line), len(content)):
                
                # get organism name for each blastx entry
                if '[' in content[i]:
                    if ']' in content[i]:
                        organism = content[i].split('[')[1][:-2]
                    else:
                        organism = content[i].split('[')[1][:-2] + ' ' + content[i+1][:-2]
                    
                    # get the uniprot acession number for the first entry of the specified species
                    if organism == species:
                        for j in range(i, -len(content), -1):
                            if 'RecName:' in content[j] and '.' in content[j]:
                                uniprot_acession = content[j].split('.')[0][1:]
                                break

                        # get the e-value and identity of the blastx mapping for the specified species entry
                        values = getBlastxValues(content, i)
                        
                        # if entry e-value and identity passes threshold, save them in a dictionary with the uniprot acession number as key
                        if float(values['evalue']) <= 0.001 and int(values['identities']) >= 40:
                            uniprot_hits[uniprot_acession] = {'evalue': values['evalue'], 'identities': values['identities']}
                            break
            
        if not line:
            # end of file - this should not happen
            print("Could not find the description section")
            return None

    return uniprot_hits


def getUniprotAnnotation(uniprot_hit, species_prefix):

    """
    Searches Uniprot for the annotation of a specific acession number

    Parameter:
        uniprot acession number
        prefix of the species respective to the acession number (e.g. _ARATH)

    Return: dictionary with the annotation of a specific acession number (TAIR, gene ID, protein name, KEGG and GO terms)
    """

    protein = ''
    gene = ''
    tair = ''
    kegg = []
    goTerms = []
    
    try:
        
        # get the Uniprot data for the the specified acession number
        htmlResponse_info = mf.urlOpenWithRetry('https://www.uniprot.org/uniprot/' + uniprot_hit + '.txt')
        stringLines_info = htmlResponse_info.data.decode('utf-8').split('\n')

        if stringLines_info != []:
            for line in stringLines_info:
                
                # check if the uniprot entry corresponds to the specified species
                if line.startswith('ID   ') and species_prefix not in line:
                    return None
                
                # get the protein name
                if line.startswith('DE   RecName: Full='):
                    protein = line.split('Full=')[1].split('{')[0].strip(' ;\n')
                
                # get the gene ID
                if gene == '':
                    if line.startswith('GN'):
                        gene = line.split('=')[1].split('{')[0].split(';')[0]
                
                # get the KEGG annotations
                if line.startswith('DR   KEGG'):
                    kegg.append(line.split(';')[1].strip(' '))
                
                # get TAIR (if species is Arabidopsis thaliana)
                if line.startswith('DR   TAIR'):
                    tair = line.split(';')[-1][1:-2]
                
                # get the GO annotations
                if line.startswith('DR   GO'):
                    line = line.split(';')
                    goTerms.append(line[2][1:3] + line[1].strip(' '))
    
    except:
        return
                
    return {'tair': tair, 'gene': gene, 'protein': protein, 'kegg': (';').join(kegg), 'goTerms': (';').join(goTerms)}


def writeDegsAnnotation(allDegs, blastx_dir, filename_out):

    """
    Write gene Uniprot annotations in a tab-delimited text file (.tab)

    Parameter:
        set with all DEGs from all comparisons (DE)
        directory with the .blastx files
        filename to output

    Output: .tab file with genes Uniprot annotations
    """
    
    degsAnnotation = {}
    
    # write file header with genes attributes
    with open(filename_out, 'w') as outfile:
        outfile.write('trinity ID\ttair\t\tprotein\tkegg\tgoTerms\te-value\tidentities\n')

        # get DEGs Uniprot acession number from blastx
        for gene in allDegs:
            filename = blastx_dir + gene + '_blastx'
            uniprot_hits = getUniprotHits(filename)

            if uniprot_hits == {}:
                print(filename, 'doesnt exist!')
                continue
            
            # get DEGs Uniprot annotations
            if uniprot_hits != None:
                for hit in uniprot_hits:
                    degAnnotation = getUniprotAnnotation(hit)
                    
                    if degAnnotation != None:
                        # add blastx stats results
                        degAnnotation['evalue'] = uniprot_hits[hit]['evalue']
                        degAnnotation['identities'] = uniprot_hits[hit]['identities']
                        degsAnnotation[gene] = degAnnotation

                        # write DEGs Uniprot annotations and blastx stats results
                        outfile.write(gene)
                        for atribute in degsAnnotation[gene]:
                            outfile.write('\t' + degsAnnotation[gene][atribute])
                        outfile.write('\n')
                        break

                    else:
                        continue

                if degAnnotation == None:
                    print(gene, 'not annotated!')
                    continue


writeDegsAnnotation(mf.getAllDegsSet(sys.argv[1]), sys.argv[2], sys.argv[3])