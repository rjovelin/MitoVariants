# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 18:31:43 2016

@author: Richard
"""

# use this script to compare mitoseek outputs for recalibrated and non-calibrated bams

# place the script in the parent folder of the subfolders with mitoseek outputs

import os


# create dicts {individual: set(snp position)}
nocal, recal = {}, {}

folders = ['./no_calibration/', './with_calibration/']
# loop over folders
for i in range(len(folders)):
    # create a list of mitoseek subfolders
    subfolders = [j for j in os.listdir(folders[i]) if 'GRCH37' in j]
    print(len(subfolders))
    # loop over each subfolder
    for k in subfolders:
        # open heteroplasmy file
        if 'mito1_heteroplasmy.txt' in os.listdir(folders[i] + k):
            # get participant
            participant = k[:k.index('_')]
            print(participant)
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
                    snp = int(line[1]) - 1
                    # check which dict to populate        
                    if i == 0:
                        # initialize dict with set
                        if participant not in nocal:
                            nocal[participant] = set()
                        else:
                            nocal[participant].add(snp)
                    elif i == 1:
                        # initialize dict with set
                        if participant not in recal:
                            recal[participant] = set()
                        else:
                            recal[participant].add(snp)
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
    print(i, 'recal:', len(recal[i]), 'nocal', len(nocal[i]), len(recal[i]) >= len(nocal[i]), len(recal[i].intersection(nocal[i])))

    
