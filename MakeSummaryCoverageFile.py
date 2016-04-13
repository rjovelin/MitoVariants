# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 18:12:58 2016

@author: RJovelin
"""


# use this script to generate a summary file with average read depth per individual

# place this script in the folder with the mitoseek subfolders

# usage: python3 MakeSummaryCoverageFile.py [options]
# - cancer_name
# - [WGS/RNAseq]: sample type
# [Refs/GRCH37]: if folder contains individuals with multiple refs or only GRCH37
# - outputfile

import os
import sys
import numpy as np
from mito_mutations import *

# get cancer name
cancer_name = sys.argv[1]
# get sample type
sample_type = sys.argv[2]
# get ref for individuals in folder
reference = sys.argv[3]
# get outputfile name
outputfile = sys.argv[4]


# check if file already exists
try:
    newfile = open(outputfile, 'r')
except:
    # write header if file does not exist
    newfile = open(outputfile, 'w')
    newfile.write('\t'.join(['participant', 'mean_read_depth', 'median_read_depth', 'tumor', 'sample_type']) + '\n')
    print('content written to new file')
else:
    # if it exists, append content
    newfile = open(outputfile, 'a')
    print('content appended to existing file')


# check if folders contains individuals with different refs or only GRCH37
if reference == 'Refs':
    # need to filter individuals not mapped to GRCH37
    subfolders = GetValidReference('./')
elif reference == 'GRCH37':
    # only indivuals mapped to GRCH37 are included in current folder
    subfolders = [i for i in os.listdir() if 'GRCH37' in i or 'RNASEQ' in i]
    # removes files if files are inluded
    to_delete = []
    for i in subfolders:
        try:
            os.listdir(i)
        except:
            to_delete.append(i)
    if len(to_delete) != 0:
        for i in to_delete:
            subfolders.remove(i)

print('# subfolders', len(subfolders))            
# loop over subfolders
for i in subfolders:
    # get the participant ID
    participant = i[:i.index('_')]
    # get the position-read depth for that participant
    reads = GetReadDepth(i + '/mito1_basecall.txt')
    # get the read counts
    ReadCounts = [reads[j] for j in reads]
    # compute mean and median read counts
    mean_reads = np.mean(ReadCounts)
    median_reads = np.median(ReadCounts)
    # write content to file
    newfile.write('\t'.join([participant, str(mean_reads), str(median_reads), cancer_name, sample_type]) + '\n')

# close file after writing
newfile.close()            
            
    












