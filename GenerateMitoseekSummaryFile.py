# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:28:27 2016

@author: RJovelin
"""

# use this script to generate a summary file with heteroplasmic variants per tumor sample
# place this script in the folder containing the mitoseek subdirectories

# usage python3 GenerateMitoseekSummaryFile.py options

# [cancer]: name of cancer (even for normal tissue)
# [normal/tumor]: tissue type
# [WGS/RNA]: use RNAseq or WGS mitoseek outputs
# BlacklistFile: file with blacklisted individuals
# MinimumReadDepth: keep individuals with median read depth > MinimumReadDepth
# PositionReadDepth: keep positions with read depth > PositionReadDepth
# outputfile

import os
import sys
import numpy as np
from mito_mutations import *

cancer = sys.argv[1]
tissue = sys.argv[2]
dataset = sys.argv[3]
BlacklistFile = sys.argv[4]
MinimumReadDepth = int(sys.argv[5])
PositionReadDepth = int(sys.argv[6])

# build outputfile name
outputfile = 'HeteroplasmySummary_' + cancer + '_' + tissue + '_' + dataset + '.txt' 


# check data sets, filter out individuals mapped to ref different than GRCH37 if WGS
if dataset == 'WGS':
    # get the valid set of subfolders using only individuals with GRCH37 reference
    folders = GetValidReference('./')
elif dataset == 'RNA':
    # bams were already selected with ref GRCH37
    folders = os.listdir('./')
    to_delete = set()
    for i in folders:
        try:
            os.listdir(i)
        except:
            to_delete.add(i)
        if 'RNASEQ' not in i:
            to_delete.add(i)
    for i in to_delete:
        folders.remove(i)
print(len(folders))        


# create a dict {participant: median read depth}
ReadDepth = MedianReadDepthIndividual(folders)

# get the set of blaclisted individuals
blacklisted = BlackListed(BlacklistFile)

# get the sample size for each position {position: set(individuals)}
# this function already ignores duplicate individuals
# and counts individuals with minimum read depth at given position
sample_size = SampleSize('./', tissue, dataset, PositionReadDepth)
print('positins with sample size', len(sample_size))
# remove blaclisted individuals and individuals with too low read depth
for position in sample_size:
    to_delete = set()
    for participant in sample_size[position]:
        if participant in blacklisted or ReadDepth[participant] <= MinimumReadDepth:
            to_delete.add(participant)
    for participant in to_delete:
        sample_size[position].remove(participant)
# delete positon with no individuals
to_delete = []
for position in sample_size:
    if sample_size[position] == 0:
        to_delete.append(position)
if len(to_delete)  != 0:
    for position in to_delete:
        del sample_size[position]
print('positions with sample size after removing low coverage individuals', len(sample_size))

# get the gene annotations
MT_annotation =  MitoAnnotation('rCRS_genes_MT.text.txt')
print('MT regions', len(MT_annotation))

# create a dictionary to store the variant info for all individuals
MitoVariants = {}

# loop over each subfolder, get the variant info, add to the MitoVariants dict    
for subfolder in folders:
    print(subfolder)
    # check that heteroplasmy file is in subfolder
    if 'mito1_heteroplasmy.txt' in os.listdir(subfolder):
        # mito1 is always the file for the tumor whether mitoseek was run with paired or individual samples
        variants = GetIndividualTumorHeteroplasmies(subfolder + '/mito1_heteroplasmy.txt', sample_size, MT_annotation)
        # update MitoVariants with variant info
        for position in variants:
            # compute the read coverage for the given position
            reads = sum(list(map(lambda x: int(x), variants[position][11:19])))
            assert type(reads) == int, 'sum of read counts should be an integer'
            # do not consider blacklisted participant IDs and participants with median read depth <= minimum
            # do not consider positions with read depth <= minimum
            if variants[position][0] not in blacklisted and ReadDepth[variants[position][0]] > MinimumReadDepth and reads > PositionReadDepth:
                if position in MitoVariants:
                    MitoVariants[position].append(variants[position])
                else:
                    MitoVariants[position] = [variants[position]]

# open file for writing
newfile = open(outputfile, 'w')

# create a list of positions in MitoVariants
positions = [i for i in MitoVariants]
positions.sort()

# write header to file
header = ['Position', 'Participant', 'Gene', 'Orientation', 'Exon_effect',
          'AA_change', 'Reference', 'Major', 'Minor', 'Major_count', 'Minor_count',
          'Sample_size', 'Forward_A', 'Forward_T', 'Forward_C', 'Forward_G', 
          'Reverse_A', 'REverse_T', 'Reverse_C', 'Reverse_G', 'FisherPval']
newfile.write('\t'.join(header) + '\n')

# loop over sorted positions
for i in positions:
    # write info for each individual at the same position to file
    for j in MitoVariants[i]:
        # convert positions back to 1-based
        newfile.write(str(i+1) + '\t' + '\t'.join(j) + '\n')

newfile.close()        
        
        
        
        
        
        