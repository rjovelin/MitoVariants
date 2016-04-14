# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 15:12:35 2016

@author: RJovelin
"""

# use this script to generate a summary file with sample frequencies of multiallelic sites
# for each tumor type

# usage: python3 MakeSummaryMultiAllelicFile.py [options]
# - HeteroplasmyFile: summary file with heteroplasmies for a given sample
# - threshold: percent of reads supporting alternative alleles
# - outputfile

import os
import sys
from mito_mutations import *

# get heteroplasmy file
HeteroplasmyFile = sys.argv[1]
# get threshold
threshold = float(sys.argv[2])
# get outputfile name
outputfile = sys.argv[3]

# extract the tumor name from the Heteroplasmy file name
cancer = HeteroplasmyFile[HeteroplasmyFile.index('_') + 1 : HeteroplasmyFile.index('_', HeteroplasmyFile.index('_') + 1)] 

# check if file already exists
try:
    newfile = open(outputfile, 'r')
except:
    # write header if file does not exist
    newfile = open(outputfile, 'w')
    newfile.write('\t'.join(['tumor', 'individual', 'position', 'gene', 'N_alleles', 'sample_size', 'sample_frequency']) + '\n')
    print('content written to new file')
else:
    # if it exists, append content
    newfile = open(outputfile, 'a')
    print('content appended to existing file')


# identify multiallelic SNPs
# get the number of alleles at each position for each individual 
# {individual : {position : [allele1, allele2, allele3, allele4]}}
AlleleCounts = IdentifyMultiAllelicPositions(HeteroplasmyFile, threshold)

# create a dict {position: [gene, sample_size]}
annotation = {}
infile = open(HeteroplasmyFile, 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.spplit('\t')
        # get position 0-based
        if int(line[0]) -1 not in annotation:
            annotation[int(line[0]) -1] = [line[2], int(line[11])]
infile.close()

# loop over individuals
for individual in AlleleCounts:
    # loop over positions
    for position in AlleleCounts[individual]:
        # only consider multiallelic snps (> 2 alleles)
        if len(AlleleCounts[individual]) > 2:
            # 
        
# report the position, alleles detected, N_alleles, gene, position_in_gene, sample_freq
    
    
    # include individual and then parse file to make another summary table?
    # include only counts, and then parse file to identify positions with multiple allelic
    # combinations and and different number of alleles per individuals?
    
    
    
newfile.close()