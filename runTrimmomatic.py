#!/usr/bin/env python3
# coding: utf-8

import os
import re


'''Runs Trimmomatic either with single-end (SE) or paired-end (PE) reads'''


def trimmomatic_PE_reads():

    """
    Runs Trimmomatic with paired-end (PE) reads, processing them according to the parameters set by the user

    Output: two .fastq.gz files (paired and unpaired reads), one for each read
    """

    while len(f_list) >= 2:
        # saves reads prefixes for each pair and sample name
        r1_pref = re.sub('fastq.gz$', '', f_list[0])
        r2_pref = re.sub('fastq.gz$', '', f_list[1])
        r_pref = r1_pref[:-3]
        
        # saves full path to each pair read and saves it in a string (space separated)
        r1_in = f_path + f_list[0]
        r2_in = f_path + f_list[1]
        readsIn = r1_in + ' ' + r2_in

        # saves full path to each output file (trimP and trimU) for each pair and saves it in a string (space separated)
        r1_trimP = outdir + r1_pref + trimP_suf
        r1_trimU = outdir + r1_pref + trimU_suf
        r2_trimP = outdir + r2_pref + trimP_suf
        r2_trimU = outdir + r2_pref + trimU_suf
        readsOut = r1_trimP + ' ' + r1_trimU + ' ' + r2_trimP + ' ' + r2_trimU
        
        # run trimmomatic for the saved input and output files
        os.system('trimmomatic ' + r_type + ' ' + readsIn + ' ' + readsOut + ' ' + parameters + ' &> ' + outdir + r_pref + '_trimmOut')
        
        # removes the already trimmed file names from the list
        f_list = f_list[2:]


def trimmomatic_SE_reads():

    """
    Runs Trimmomatic with single-end (SE) reads, processing them according to the parameters set by the user

    Output: two .fastq.gz files (paired and unpaired), one for each input read
    """

    for f in f_list:
        # saves read prefix
        r_pref = re.sub('fastq.gz$', '', f)

        # save full path to the read
        readIn = f_path + f
        
        # save full path to the output files (trimP and trimU, space separated)
        r_trimP = outdir + r_pref + trimP_suf
        r_trimU = outdir + r_pref + trimU_suf
        readsOut = r_trimP + ' ' + r_trimU
        
        # run trimmomatic for the saved input and output files
        os.system('trimmomatic ' + r_type + ' ' + readIn + ' ' + readsOut + ' ' + parameters + ' &> ' + outdir + r_pref + '_trimmOut')


#get inputs from user
r_type = str(input('reads type? (PE/SE) '))
f_path = str(input('path to directory of files to be trimmed? '))
outdir = str(input('output directory path? '))
parameters = str(input('parameters command? (see Trimmomatic manual) '))

#creates a list of all file names to be trimmed and sorts them alphabetically
f_list = os.listdir(f_path)
f_list = [file for file in f_list if file.endswith('fastq.gz')]
f_list.sort()

# define output files suffixes
trimP_suf = 'trimP.fastq.gz'
trimU_suf = 'trimU.fastq.gz'


#if the reads are pair ended, run trimmomatic for each pair of reads (1/2 or left/right)
if r_type == 'PE':
    trimmomatic_PE_reads()

#if the reads are single ended, run trimmomatic for each read
else:
    trimmomatic_SE_reads()