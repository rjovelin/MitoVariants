# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 16:03:58 2016

@author: RJovelin
"""


# use this script to check if mitoseek identifies > 2 alleles in the summary file
# redirects to stdout
# place this script in the folder containing the mitoseek subdirectories

# usage python3 CheckMultiAlleles.py 
# - [WGS/RNA]: if RNAseq or WGS data is used
# tumor_type

import os
import sys


# get tumor type
tumor = os.getcwd()
tumor = tumor[tumor.index('mitoseek_hp1_mbq13')+len('mitoseek_hp1_mbq13')+1:]
tumor = tumor[:tumor.index('/')]

folders = os.listdir('./')
to_delete = []
for i in folders:
    try:
        os.listdir(i)
    except:
        to_delete.append(i)
for i in to_delete:
    folders.remove(i)

# make a set of valid nucleotides
valid_bases = {'A', 'T', 'C', 'G'}
print('individual', 'ref', 'major', 'minor', sep = '\t')

# loop over each subfolder 
for subfolder in folders:
    # get individual ID
    individual = subfolder[:subfolder.index('_')]
    # open file for reading
    infile = open(subfolder + '/mito1_heteroplasmy.txt', 'r')
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get position 1-based
            position = line[1]
            ref, major, minor = line[2], line[14], line[15]
            # get major and minor alleles
            if major.upper() not in valid_bases or minor.upper() not in valid_bases:
                print(individual, position, ref, major, minor, sep = '\t')
    infile.close()        
    
