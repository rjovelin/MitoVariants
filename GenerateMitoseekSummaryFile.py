# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:28:27 2016

@author: RJovelin
"""

# use this script to generate a summary file with heteroplasmic variants per tumor sample
# place this script in the folder containing the mitoseek subdirectories

# usage python3 GenerateMitoseekSummaryFile.py options
# [tumor/normal]: consider normal or tumor sample
# [WGS/RNA]: use RNAseq or WGS mitoseek outputs
# outputfile

import os
import sys
from somatic_mutations import *

tissue = sys.argv[1]
dataset = sys.argv[2]
outputfile = sys.argv[3]

# get the sample size for each position {position: N}
# this function already ignores duplicate individuals
sample_size = SampleSize('./', tissue, dataset)
print(len(sample_size))

# check data sets, filter out individuals mapped to ref different than GRCH37 if WGS
if dataset == 'WGS':
    # get the valid set of subfolders using only individuals with GRCH37 reference
    folders = GetValidReference('./')
elif dataset == 'RNA':
    # bams were already selected with ref GRCH37
    folders = os.listdir('./')
    to_delete = []
    for i in folders:
        try:
            os.listdir(i)
        except:
            to_delete.append(i)
    for i in to_delete:
        folders.remove(i)
print(len(folders))        

# get the gene annotations
MT_annotation =  MitoAnnotation('rCRS_genes_MT.text.txt')
print(len(MT_annotation))

# create a dictionary to store the variant info for all individuals
MitoVariants = {}

# loop over each subfolder, get the variant info, add to the MitoVariants dict    
for subfolder in folders:
    print(subfolder)
    # mito1 is always the file for the tumor whether mitoseek was run with paired or individual samples
    variants = GetIndividualTumorHeteroplasmies(subfolder + '/mito1_heteroplasmy.txt', sample_size, MT_annotation)
    # update MitoVariants with variant info
    for position in variants:
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
        
        
        
        
        
        