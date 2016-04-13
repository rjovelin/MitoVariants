# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 13:52:57 2016

@author: RJovelin
"""


# use this function to filter out positions that are variable in RNAs and in DNA

import os
import sys
from somatic_mutations import *

# usage python Filter_RNA_Heteroplasmies.py options
# - RNA_summary_file: summary file with heteroplasmies in RNAs
# - DNA_summary_file: sumary file with heteroplasmies in DNA
# - RNA_folder: folder with RNA mitoseek outputs
# - WGS_folder: folder with WGS mitoseek outputs
# - suffix: suffix of the participant ID in the subfolder name
# - [soft/hard]: filter to compare variable positions in RNA and DNA.
#                Remove all variable positions in DNA (hard)
#                Or remove variable positions in DNA if alleles are the same in RNA and DNA
# - minimum_coverage: minimum read depth to keep positions               
# - outputfile


# get the RNA summary file
SummaryRNA = sys.argv[1]
# get the DNA summary file
SummaryDNA = sys.argv[2]
# get the folder with mitoseek RNA outputs
RNA_folder = sys.argv[3]
# get the folder with mitoseek WGS outputs
WGS_folder = sys.argv[4]
# get the suffix
suffix = sys.argv[5]
# get the strength of the RNA filter
RNAFilter = sys.argv[6]
# get the minimum read depth
minimum_coverage = int(sys.argv[7])
# get the outputfile name
outputfile = sys.argv[8]


# parse the summary files into a dict of dict {participantID: {position : [information]} 
RNA_snps = GetVariablePositions(SummaryRNA)
WGS_snps = GetVariablePositions(SummaryDNA)
print('RNA_snps', len(RNA_snps))
print('DNA_snps', len(WGS_snps))

# remove participants without RNA and WGS data
RNA_snps, WGS_snps = RemoveUniqueParticipants(RNA_snps, WGS_snps)
print('RNA_snps after removing unique participants', len(RNA_snps))

print(os.getcwd())
# move to directory containing the subfolders of the mitoseek WGS outputs
os.chdir(WGS_folder)
print(os.getcwd())

# remove positions with RNA variants that have no coverage in DNA
RNA_snps = RemovePositionWithoutCoverage(RNA_snps, suffix, minimum_coverage)
print('RNA_snps after removing positions without coverage', len(RNA_snps))

# move to the directory containing the subfolders of the mitoseek RNA outputs
os.chdir(RNA_folder)
print(os.getcwd())


# correct sample size at each position
# get the sample size of individuals with RNAseq and WGS that have minimum coverage in RNA at each position
# make a list of subfolders from RNA mitoseek outputs
RNASubFolders = [i for i in os.listdir(RNA_folder)]
to_delete = []
for i in RNASubFolders:
    try:
        os.listdir(i)
    except:
        to_delete.append(i)
    else:
        if i[:i.index('_TP_RNASEQ')] not in RNA_snps:
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
for ID in RNA_snps:
    # loop over position
    for position in RNA_snps[ID]:
        # get the sample size
        RNA_snps[ID][position][11] = str(sample[position])


# filter variants observed in DNA
# loop over participants in RNA
for participant in RNA_snps:
    # create a list of positions to remove
    to_delete = []
    # loop over positions
    for position in RNA_snps[participant]:
        # check if position is recorded in WGS
        if position in WGS_snps[participant]:
            # check filter strength
            if RNAFilter == 'soft':
                # remove positions for which variants are the same in RNA and DNA
                # keep variable positions in DNA if alleles are different in RNA and DNA
                rna_alleles = set(list(map(lambda x: x.upper(), RNA_snps[participant][position][7:9])))
                dna_alleles = set(list(map(lambda x: x.upper(), WGS_snps[participant][position][7:9])))
                if rna_alleles == dna_alleles:
                    to_delete.append(position)
            elif RNAFilter == 'hard':
                # remove position
                to_delete.append(position)
    if len(to_delete) != 0:
        for position in to_delete:
           del RNA_snps[participant][position]

# remove participants with no variants
to_delete = []
for i in RNA_snps:
    if len(RNA_snps[i]) == 0:
        to_delete.append(i)
for i in to_delete:
    del RNA_snps[i]
print('RNA_snps after filtering variants in DNA', len(RNA_snps))


# create a dict with position as key and list of lists with info for all participants
# {position : [[information_ID1], [information_ID2]}
RNA_variants = {}
for ID in RNA_snps:
    # get positions
    for position in RNA_snps[ID]:
        if position in RNA_variants:
            RNA_variants[position].append(RNA_snps[ID][position])
        else:
            RNA_variants[position] = [RNA_snps[ID][position]]
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

