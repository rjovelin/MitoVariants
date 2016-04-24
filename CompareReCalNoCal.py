# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 18:31:43 2016

@author: Richard
"""

# use this script to compare mitoseek outputs for recalibrated and non-calibrated bams

# place the script in the parent folder of the subfolders with mitoseek outputs

import os


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
                    
print('recal', len(recal))
print('nocal', len(nocal))

# remove unique participant
to_delete = [i for i in recal if i not in nocal]
print('remove from recal', len(to_delete))
if len(to_delete) != 0:
    for i in to_delete:
        del recal[i]
to_delete = [i for i in nocal if i not in recal]
print('remove from nocal', len(to_delete))
if len(to_delete) != 0:
    for i in to_delete:
        del nocal[i]
    
print('recal', len(recal))
print('nocal', len(nocal))

# compare snps in same individuals
for i in recal:
    recal_snps = set(j for j in recal[i])
    nocal_snps = set(j for j in nocal[i])
    print(i, 'recal:', len(recal[i]), 'nocal', len(nocal[i]), len(recal[i]) >= len(nocal[i]), len(recal_snps.intersection(nocal_snps)))

print('\n')
print('======\n' * 2)


for i in recal: 
    for j in recal[i]:
        if j in nocal[i]:
            print(i, recal[i][j][0] == nocal[i][j][0], recal[i][j][3] == nocal[i][j][3], recal[i][j][2], nocal[i][j][2], recal[i][j][5], nocal[i][j][5])
           
