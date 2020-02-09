#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

#gets inputs from user
f_path = str(input('Path to directory with files to analyse?'))
outdir = str(input('Output directory?'))
#create output directory if it does not exist
if not os.path.isdir(outdir):
	os.system('mkdir ' + outdir)

def fastq_analysis(tool):

	if tool == 'multiqc':

		#run multiqc for all files of a given directory
		os.system('multiqc ' + f_path + ' -o ' + outdir)

	else:

		#save filenames from a given directory to analyse in a list
		f_list = os.listdir(f_path)
		f_list.sort()


		if tool == 'fastqc':

			#run fastqc for each file of the given directory
			for f in f_list:
				os.system('fastqc ' + f_path + f + ' -o ' + outdir)


		if tool == 'fastqscreen':

			#get inputs from user
			aligner = str(input('Aligner? (BWA/Bowtie/Bowtie2)'))
			built_indices = str(input('Download pre-build aligner indices? (y/n)'))
			fqs_path = str(input('Path to fastqc screen directory? (excl. exec)'))
			parameters = str(input('Fastq screen parameters command? (see Fastq screen manual)'))

			if built_indices == 'y':
				gen_outdir = str(input('Genomes output directory?'))
				#pre-built aligner indices of commonly used genomes, downloaded directly from the Babraham Bioinformatics website
				os.system('fastq_screen --get_genomes -o ' + gen_outdir)
				#copy the new configuration file to fastqscreen directory
				os.system('cp gen_outdir/fastq_screen.conf ' + fqs_path)

			#tag and/or filter reads that don't align with the searched reference genomes
			for f in f_list:
				os.system(fqs_path + 'fastq_screen --aligner ' + aligner + ' ' + parameters + ' ' + f_path + f)

tools = str(input('Which tools? (fastqc,fastqscreen,multiqc - comma separated)'))
tools = tools.split(',')

for tool in tools:
	fastq_analysis(tool)