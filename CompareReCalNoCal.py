# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 18:31:43 2016

@author: Richard
"""

# use this script to compare mitoseek outputs for recalibrated and non-calibrated bams

# place the script in the parent folder of the subfolders with mitoseek outputs

import os
import numpy as np

# create dicts {individual: {position: [major, major_count, major_freq, minor, minor_count, minor_freq]}}
nocal, recal = {}, {}

folders = ['./no_calibration/', './with_calibration/']
# loop over folders
for i in range(len(folders)):
    # create a list of mitoseek subfolders
    subfolders = [j for j in os.listdir(folders[i]) if 'GRCH37' in j]
    # loop over each subfolder
    for k in subfolders:
        # open heteroplasmy file
        if 'mito1_heteroplasmy.txt' in os.listdir(folders[i] + k):
            # get participant
            participant = k[:k.index('_')]
            filename = folders[i] + k + '/mito1_heteroplasmy.txt'
            infile = open(filename)
            # skip header
            infile.readline()
            # loop over file
            for line in infile:
                line = line.rstrip()
                if line != '':
                    line = line.split('\t')
                    # get position index 0-based
                    position = int(line[1]) - 1
                    # get major and minor allele
                    major, minor = line[14], line[15]
                    # get read counts
                    reads = {'A': int(line[3]) + int(line[7]),
                             'T': int(line[4]) + int(line[8]),
                             'C': int(line[5]) + int(line[9]),
                             'G': int(line[6]) + int(line[10])}
                    # get major and minor counts
                    major_count, minor_count = int(line[16]), int(line[17])
                    # get major and minor frequencies
                    total = sum([reads[base] for base in reads])
                    major_freq, minor_freq = round(major_count/total, 3), round(minor_count/total, 3)
                    # check which dict to populate 
                    if i == 0:
                        # initialize dict with dict
                        if participant not in nocal:
                            nocal[participant] = {}
                        # populate dict {individual: {position: [major, major_count, major_freq, minor, minor_count, minor_freq]}}
                        nocal[participant][position] = [major, major_count, major_freq, minor, minor_count, minor_freq]    
                    elif i == 1:
                        # initialize dict with dict
                        if participant not in recal:
                            recal[participant] = {}
                        recal[participant][position] = [major, major_count, major_freq, minor, minor_count, minor_freq]
            # close file after reading
            infile.close()
                    

# remove unique participant
to_delete = [i for i in recal if i not in nocal]
if len(to_delete) != 0:
    for i in to_delete:
        del recal[i]
to_delete = [i for i in nocal if i not in recal]
if len(to_delete) != 0:
    for i in to_delete:
        del nocal[i]
    
print('\n')

print('participant', 'recal', 'N', 'nocal', 'N', 'N_recal >= N_nocal', 'N_common', sep = '\t')
# compare snps in same individuals
for i in recal:
    recal_snps = set(j for j in recal[i])
    nocal_snps = set(j for j in nocal[i])
    print(i, 'recal', len(recal[i]), 'nocal', len(nocal[i]), len(recal[i]) >= len(nocal[i]), len(recal_snps.intersection(nocal_snps)), sep = '\t')

print('\n')
print('======\n' * 2)

same_major = 0
same_minor = 0
total = 0

for i in recal: 
    for j in recal[i]:
        if j in nocal[i]:
            if recal[i][j][0] == nocal[i][j][0]:
                same_major += 1
            if recal[i][j][3] == nocal[i][j][3]:
                same_minor += 1
            total += 1

print('allele', 'N_same_alleles', 'freq', sep = '\t')
print('major', same_major, round(same_major / total, 3), sep = '\t')
print('minor', same_minor, round(same_minor / total, 3), sep = '\t')            


print('\n')
print('======\n' * 2)


# print median, mean for the difference in frequencies for major and minor alleles 
MajorFreqDiff, MinorFreqDiff = [], []
for i in recal: 
    for j in recal[i]:
        if j in nocal[i]:
            MajorFreqDiff.append(abs(round(recal[i][j][2] - nocal[i][j][2], 6)))
            MinorFreqDiff.append(abs(round(recal[i][j][5] - nocal[i][j][5], 6)))

print('allele', 'median', 'mean', sep = '\t')
print('major diff', round(np.median(MajorFreqDiff), 4), round(np.mean(MajorFreqDiff), 4), sep = '\t')
print('minor diff', round(np.median(MinorFreqDiff), 4), round(np.median(MinorFreqDiff), 4), sep = '\t')

print('\n')
print('======\n' * 2)
         
### check read depth of unique snps

# create dicts {participant: set(unique snps)}
recal_unique, nocal_unique = {}, {}

for i in recal:
    for j in recal[i]:
        if j not in nocal[i]:
            if i not in recal_unique:
                recal_unique[i] = set()
            recal_unique[i].add(j)

for i in nocal:
    for j in nocal[i]:
        if j not in recal[i]:
            if i not in nocal_unique:
                nocal_unique[i] = set()
            nocal_unique[i].add(j)

# create lists to store the total read depth for each unique position
TotalReadDepthRecal, TotalReadDepthNocal = [], []
# create a dict: {participant: {position: read depth}}
ReadDepthRecal, ReadDepthNocal = {}, {}
# loop over participant in recal_unique
for i in recal_unique:
    # open file for reading
    infile = open('./with_calibration/' + i + '_recal_TP_GRCH37/mito1_basecall.txt')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            pos = int(line[1]) -1
            reads = sum(list(map(lambda x:int(x), line[3:])))
            if i not in ReadDepthRecal:
                ReadDepthRecal[i] = {}
            ReadDepthRecal[i][pos] = reads
    infile.close()

for i in nocal_unique:
    # open file for reading
    infile = open('./no_calibration/' + i + '_nocal_TP_GRCH37/mito1_basecall.txt')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            pos = int(line[1]) -1
            reads = sum(list(map(lambda x:int(x), line[3:])))
            if i not in ReadDepthNocal:
                ReadDepthNocal[i] = {}
            ReadDepthNocal[i][pos] = reads
    infile.close()
    

for i in recal_unique:
    for j in recal_unique[i]:
        TotalReadDepthRecal.append(ReadDepthRecal[i][j])

for i in nocal_unique:
    for j in nocal_unique[i]:
        TotalReadDepthNocal.append(ReadDepthNocal[i][j])

print('calibration', 'N', 'median', 'mean', 'min', 'max', sep = '\t')        
print('recal', len(TotalReadDepthRecal), round(np.median(TotalReadDepthRecal), 3), round(np.mean(TotalReadDepthRecal), 3), min(TotalReadDepthRecal), max(TotalReadDepthRecal), sep = '\t')
print('nocal', len(TotalReadDepthNocal), round(np.median(TotalReadDepthNocal), 3), round(np.mean(TotalReadDepthNocal), 3), min(TotalReadDepthNocal), max(TotalReadDepthNocal), sep = '\t')

print('\n')
print('======\n' * 2)

TotalReadDepthNocal.sort()
TotalReadDepthRecal.sort()
# count number of position with read depth higher and below potential cutoff (1000 reads)
higher_recal, lower_recal, higher_nocal, lower_nocal = 0, 0, 0, 0
for i in TotalReadDepthRecal:
    if i < 1000:
        lower_recal += 1
    elif i >= 1000:
        higher_recal += 1

for i in TotalReadDepthNocal:
    if i < 1000:
        lower_nocal += 1
    elif i >= 1000:
        higher_nocal += 1
    
print('calibration', 'depth', 'N', 'freq', sep = '\t')
print('recal', 'lower than 1000', lower_recal, round(lower_recal / len(TotalReadDepthRecal), 3) * 100, sep = '\t')
print('recal', 'higher than 1000', higher_recal, round(higher_recal / len(TotalReadDepthRecal), 3) * 100, sep = '\t')
print('nocal', 'lower than 1000', lower_nocal, round(lower_nocal / len(TotalReadDepthNocal), 3) * 100, sep = '\t')
print('nocal', 'higher than 1000', higher_nocal, round(higher_nocal / len(TotalReadDepthNocal), 3) * 100, sep = '\t')



