# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 10:49:31 2016

@author: RJovelin
"""

# use this script to save a dictionary {participant: [[read_depth_snps], [read_depth_not_variable]]}

import os
import sys
import json
from mito_mutations import *

# place this script in the folder with mitoseek subfolders output

# usage GetReadDepthMutations.py [options]
# - HeteroSummaryFile: file with heteroplasmies
# -[Refs/GRCH37]: if folder contains individuals with multiple refs or only GRCH37
# - outputfile: file to store the data setructure in json format


# get heteroplasmy file
HeteroSummaryFile = sys.argv[1]
# get reference
reference = sys.argv[2]
# get the outputfile name
outputfile = sys.argv[3] 

# get the positions with heteroplasmies for each individual
# create a dict {participant: positions}
VariablePositions = {}
# open file for reading
infile = open(HeteroSummaryFile, 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get participant ID
        participant = line[1]
        # get position 0-based
        position = int(line[0]) -1
        # populate dict
        if participant not in VariablePositions:
            VariablePositions[participant] = set()
        VariablePositions[participant].add(position)
infile.close()


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


# create dict {participant: [[read_depth_snps], [read_depth_not_variable]]}
ReadDepth = {}
# loop over subfolders
for i in subfolders:
    # get the participant ID
    participant = i[:i.index('_')]
    # get the position-read dpeth for that participant
    reads = GetReadDepth(i + '/mito1_basecall.txt')
    # initialize dict 
    if participant not in ReadDepth:
        ReadDepth[participant] = [[], []]
    # check if position is variable or not
    for position in reads:
        if position in VariablePositions[participant]:
            # position is variable, append read depth to first list
            ReadDepth[participant][0].append(reads[position])
        else:
            ReadDepth[participant][1].append(reads[position])
            
# save dictionary structure as a json file
newfile = open(outputfile, 'w')
# dump dict to file
json.dump(ReadDepth, newfile, indent = 4, sort_keys = True)
newfile.close()


	
   
     
            


