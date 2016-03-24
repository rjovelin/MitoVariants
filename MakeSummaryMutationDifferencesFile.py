# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:21:17 2016

@author: RJovelin
"""


# use this script to make a summary file with differences between the number
# of mutations in RNA and in DNA for each individual and tumor type

# place this script in the folder with the mitoseek subfolders

# usage: python3 MakeSummaryMutationDifferencesFile.py [options]
# - cancer_name: name of tumor
# - HeteroplasmyFileRNA: summary file of mitoseek output for RNA
# - HeteroplasmyFileWGS: summary file of mitoseek output for WGS
# - outputfile

import os
import sys


# get cancer name
cancer_name = sys.argv[1]
# get heteroplasmy RNA file
HeteroplasmyRNA = sys.argv[2]
# get heteroplasmy WGS file
HeteroplasmyWGS = sys.argv[3]
# get outputfile name
outputfile = sys.argv[4]


# check if file already exists
try:
    newfile = open(outputfile, 'r')
except:
    # write header if file does not exist
    newfile = open(outputfile, 'w')
    newfile.write('\t'.join(['participant', 'Mutational_difference', 'tumor']) + '\n')
    print('content written to new file')
else:
    # if it exists, append content
    newfile = open(outputfile, 'a')
    print('content appended to existing file')

# create a dict {participant: set(variable positions)} for mutations in RNA
MutationsRNA = {}
# open heteroSummary file for reading
infile = open(HeteroplasmyRNA, 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        participant = line[1]
        position = int(line[0]) - 1
        if participant in MutationsRNA:
            MutationsRNA[participant].add(position)
        else:
            MutationsRNA[participant] = set()
            MutationsRNA[participant].add(position)
infile.close()

# create a dict {participant: set(variable positions)} for mutations in WGS
MutationsWGS = {}
# open heteroSummary file for reading
infile = open(HeteroplasmyWGS, 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        participant = line[1]
        position = int(line[0]) - 1
        if participant in MutationsWGS:
            MutationsWGS[participant].add(position)
        else:
            MutationsWGS[participant] = set()
            MutationsWGS[participant].add(position)
infile.close()

# create a dict {participant: N_mutations_RNA - N_mutations_WGS}
Variants = {}
# loop over participant in MutationsRNA
for participant in MutationsRNA:
    # consider only individuals with WGS and RNAseq data
    if participant in MutationsWGS:
        Variants[participant] = len(MutationsRNA[participant]) - len(MutationsWGS[participant])

# loop over participants having both RNA and WGS data
for participant in Variants:
    # write content to file
    newfile.write('\t'.join([participant, str(Variants[participant]), cancer_name]) + '\n')
    
# close file after writing
newfile.close()            
