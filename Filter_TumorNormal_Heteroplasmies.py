# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:10:48 2016

@author: RJovelin
"""

# use this function to filter out positions that are variable in tumor and in normal sample
# Precondition: sites have been filtered between RNA and DNA and only RNA-specific sites are reported in summary files

import os
import sys
from mito_mutations import *

# usage python Filter_TumorNormal_Heteroplasmies.py options
# - Tumor_RNAspecific_summary_file: summary file with RNA-specific heteroplasmies in tumor
# - Normal_RNAspecific_summary_file: sumary file with RNA-specific heteroplasmies in normal




# - RNA_folder: folder with RNA mitoseek outputs
# - WGS_folder: folder with WGS mitoseek outputs
# - suffix: suffix of the participant ID in the subfolder name
# - minimum_coverage: minimum read depth to keep positions               
# - outputfile


# get the Tumor summary file
SummaryTumor = sys.argv[1]
# get the Normal summary file
SummaryNormal = sys.argv[2]
# get the folder with mitoseek RNA outputs for normal sample
# get the minimum read depth
minimum_coverage = int(sys.argv[3])
# get the outputfile name



# parse the summary files into a dict of dict {participantID: {position : [information]} 
# remove participants without tumor and normal data
# remove positions with RNA-tumor variants that have no coverage (or less than threshold) in normal-tumor
# correct sample size at each position

cancer_name = SummaryTumor[SummaryTumor.index('Summary_') + len('Summary_') : SummaryTumor.index('tumor') -1]
assert cancer_name == SummaryNormal[SummaryNormal.index('Summary_') + len('Summary_') : SummaryNormal.index('normal') -1], 'cancer names do not match'

# verify that arguments are passed appropriately
assert 'tumor' in SummaryTumor and 'normal' in SummaryNormal, 'summary files should be tumor and normal'
print('QCed files matching')

# build outputfile with comand option arguments
outputfile = 'HeteroplasmySummary_' + cancer_name + '_TumorScpecific.txt'

# parse the summary files into a dict of dict {participantID: {position : [information]} 
Tumor_snps = GetVariablePositions(SummaryTumor)
Normal_snps = GetVariablePositions(SummaryNormal)
print('Tumor_snps', len(Tumor_snps))
print('Normal_snps', len(Normal_snps))

# remove participants not in Tumor and Normal
Tumor_snps, Normal_snps = RemoveUniqueParticipants(Tumor_snps, Normal_snps)
print('N individuals in Tumor after removing unique participants', len(Tumor_snps))

print(os.getcwd())
# move to directory containing the subfolders of the mitoseek RNA Normal outputs
os.chdir('NT_RNASeq/')
print(os.getcwd())

# remove positions with RNA variants that have not enough coverage in Normal RNA
Tumor_snps = RemovePositionWithoutCoverage(Tumor_snps, '_NT_RNASEQ', minimum_coverage)
print('N individuals in Tumor after removing positions without sufficient coverage', len(Tumor_snps))

# move to the directory containing the subfolders of the mitoseek RNA outputs
os.chdir('../TP_RNASeq/')
print(os.getcwd())

# correct sample size at each position
# get the sample size of individuals with RNAseq and WGS that have minimum coverage in RNA at each position
# make a list of subfolders from RNA mitoseek outputs
RNASubFolders = [i for i in os.listdir('./')]
to_delete = []
for i in RNASubFolders:
    try:
        os.listdir(i)
    except:
        to_delete.append(i)
    else:
        if i[:i.index('_TP_RNASEQ')] not in Tumor_snps:
            to_delete.append(i)
for i in to_delete:
    RNASubFolders.remove(i)
# create a dict with position 0-based {position: number of indidivuals}    
sample = {}
# loop over directories
for subfolder in RNASubFolders:
    # get basecall file
    infile = open(subfolder + '/' + 'mito1_basecall.txt', 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get position 0-based
            position = int(line[1]) - 1
            # compute coverage
            reads = sum(list(map(lambda x: int(x), line[3:])))
            assert type(reads) == int, "reads should be an integer"
            # do not consider positions with read depth < minimum
            if reads > minimum_coverage:
                # update dict with count
                sample[position] = sample.get(position, 0) + 1
    infile.close()
# loop over participant
for ID in Tumor_snps:
    # loop over position
    for position in Tumor_snps[ID]:
        # get the sample size
        Tumor_snps[ID][position][11] = str(sample[position])

# filter variants observed in Normal RNA from Tumor RNA variants
# loop over participants in Tumor
for participant in Tumor_snps:
    # create a list of positions to remove
    to_delete = []
    # loop over positions
    for position in Tumor_snps[participant]:
        # check if position is recorded in WGS
        if position in Normal_snps[participant]:
            # remove variable positions regardless of the variants themselves
            to_delete.append(position)
    print('delete {0} positions in {1}'.format(len(to_delete), participant))
    if len(to_delete) != 0:
        for position in to_delete:
           del Tumor_snps[participant][position]

# remove participants with no variants
to_delete = [i for i in Tumor_snps if len(Tumor_snps[i]) == 0]
print('delete {0} individuals'.format(len(to_delete)))
for i in to_delete:
    del Tumor_snps[i]
print('N individual in tumor after filtering variants in Normal', len(Tumor_snps))


# create a dict with position as key and list of lists with info for all participants
# {position : [[information_ID1], [information_ID2]}
RNA_variants = {}
for ID in Tumor_snps:
    # get positions
    for position in Tumor_snps[ID]:
        if position in RNA_variants:
            RNA_variants[position].append(Tumor_snps[ID][position])
        else:
            RNA_variants[position] = [Tumor_snps[ID][position]]
# create a list of positions
positions = [i for i in RNA_variants]
positions.sort()

# move to parent directory in which the script was launched
os.chdir('../')

# open file for writing
newfile = open(outputfile, 'w')
# write header to file
header = ['Position', 'Participant', 'Gene', 'Orientation', 'Exon_effect',
          'AA_change', 'Reference', 'Major', 'Minor', 'Major_count', 'Minor_count',
          'Sample_size', 'Forward_A', 'Forward_T', 'Forward_C', 'Forward_G', 
          'Reverse_A', 'REverse_T', 'Reverse_C', 'Reverse_G', 'FisherPval']
newfile.write('\t'.join(header) + '\n')

# loop over positions
for i in positions:
    # loop over each SNP at that position for all participants
    for j in RNA_variants[i]:
        # write info to file, has already position 1-based and participant ID
        newfile.write('\t'.join(j) + '\n')

newfile.close()

