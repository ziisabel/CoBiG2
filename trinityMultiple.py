#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

#creates a list with all cleaned paired reads
path = str(input('path to cleaned reads?'))
f_list = os.listdir(path)
f_list.sort()

rightStr = ''
leftStr = ''

#creates a string with only rigth reads files and another with only left reads files (comma separated)
for f in f_list:
    f = path + f
    if int(f.find('_1')) == -1:
        rightStr += f + ','
    else:
        leftStr += f + ','
rightStr = rightStr[:-1]
leftStr = leftStr[:-1]

# gets parameters inputs from user
SS_lib_type = str(input('SS_lib_type: (RF/FR/F/R)')) #RF
max_memory = str(input('max_memory: (G)')) #50G
seq_Type = str(input('seq_Type: (cfa/cfq/fa/fq)')) #fq
CPU = str(input('How many threads?')) #12
dir_name = str(input('Output directory name?')) #outTrinity

# run the analysis with the parameters set by the user
os.system('Trinity --seqType ' + seq_Type + ' --max_memory ' + max_memory + 'G --left ' + leftStr + ' --right ' + rightStr + ' --SS_lib_type ' + SS_lib_type + ' --CPU ' + CPU + ' --output ' + dir_name)

# performs and output basic assembly statistics
os.system('TrinityStats.pl Trinity.fasta >& stats.log')