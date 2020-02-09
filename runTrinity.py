#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

#creates a list with all trimmed paired file names
path = str(input("path to files?"))
f_list = os.listdir(path)
f_list.sort()

rightStr = ""
leftStr = ""

#creates a string with only rigth reads files and another with only left reads files (comma separated)
for f in f_list:
    f = path + f
    if int(f.find("1")) < 0:
        rightStr += f + ","
    else:
        leftStr += f + ","

#criar uma string para cada direcao com os nomes de todas as reads correspondentes
rightStr = rightStr[:-1]
leftStr = leftStr[:-1]

#RF
SS_lib_type = str(input("SS_lib_type? (RF/FR/F/R)"))
#50G
max_memory = str(input("max_memory? (G)"))
#fq
seq_Type = str(input("seq_Type? (cfa/cfq/fa/fq)"))
#12
CPU = str(input("How many threads (CPU)?"))
#outTrinity
dir_name = str(input("Output directory name?"))
#add other parameters

#os.system("Trinity --seqType " + seqType + " --max_memory " + max_memory + " --left " + leftStr + " --right " + rightStr + " --SS_lib_type " + SS_lib_type + "--CPU " + CPU + " --output " + dir_name)
os.system("echo Trinity --seqType " + seq_Type + " --max_memory " + max_memory + " --left " + leftStr + " --right " + rightStr + " --SS_lib_type " + SS_lib_type + " --CPU " + CPU + " --output " + dir_name)

os.system("TrinityStats.pl Trinity.fasta >& stats.log")