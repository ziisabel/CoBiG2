#!/usr/bin/env python3
import os

# load genome index to memory
genomeIndices = input('Genome indices directory: ')
os.system('STAR --genomeLoad LoadAndExit --genomeDir ' + genomeIndices)

# get files names list
f_path = input('FastQ directory (cleaned reads): ')
f_list = os.listdir(f_path)
f_list = [file for file in f_list if file.endswith('.fastq.gz')]
f_list.sort()

threads = input('Threads: ')
mappingsDir = input('Mappings directory: ')

# if PE get samples names list to group forward and reverse files
r_type = input('PE/SE: ')

if r_type == 'PE':
    
    while len(f_list) > 1:
        
        forward = f_path + f_list[0]
        reverse = f_path + f_list[1]
        print('Forward file path:', forward)
        print('Reverse file path:', reverse)
        
        sample = f_list[0].replace('_1.fastq.gz', '')
        os.system('mkdir ' + mappingsDir + '/' + sample + '_RUN')
        
        os.system('STAR --runThreadN ' + threads +
		' --readFilesIn ' + forward + ' ' + reverse +
		' --readFilesCommand zcat' +
		' --genomeDir ' + genomeIndices +
		' --outSAMtype BAM Unsorted' +
		' --quantMode GeneCounts' +
		' --outFilterMismatchNoverLmax 0.05' +
		' --outFileNamePrefix ' + mappingsDir + '/' + sample + '_RUN/' + sample)
        f_list = f_list[2:]

else:
    for file in f_list:
        filename = f_path + file
        print('Sample file path:', filename)
        sample = file.replace('.fastq.gz', '')
        
        os.system('mkdir' + mappingsDir + '/' + sample + '_RUN')
        os.system('STAR --runThreadN ' + threads +
		' --readFilesIn ' + f_path + sample +
		' --readFilesCommand zcat' +
		' --genomeDir ' + genomeIndices +
		' --outSAMtype BAM Unsorted' +
		' --quantMode GeneCounts' +
		' --outFilterMismatchNoverLmax 0.05' +
		'--outFileNamePrefix ' + mappingsDir + '/' + sample + '_RUN/' + sample)
                    
# unload genome index from memory
os.system('STAR --genomeLoad Remove --genomeDir ' + genomeIndices)
