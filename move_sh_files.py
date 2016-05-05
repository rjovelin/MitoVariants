# -*- coding: utf-8 -*-
"""
Created on Wed May  4 14:44:26 2016

@author: RJovelin
"""

# use this script to move sh files to corresponding folders
# place this script in mitoseek_hp1_mbq20

import os

# create a set of tumor from input file names
inputfiles = [i for i in os.listdir() if '.txt' in i and ('WGS' in i or 'RNAseq' in i)]

tumor_types = set()
for filename in inputfiles:
    tumor = filename[:filename.index('_')]
    tumor_types.add(tumor)

# create a dict {datatype : {tumor: [individuals]}}
MapTumor = {}
for filename in inputfiles:
    if 'WGS' in filename:
        datatype = 'WGS'
    elif 'RNAseq' in filename:
        datatype = 'RNAseq'
    # get tumor name
    tumor = filename[:filename.index('_')]
    # check if datatype in dict
    if datatype not in MapTumor:
        MapTumor[datatype] = {}
    # check if tumor in inner dict
    if tumor not in MapTumor[datatype]:
        MapTumor[datatype][tumor] = []
    # open file for reading
    infile = open(filename)
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            individual = line[0]
            MapTumor[datatype][tumor].append(individual)
    infile.close()

# create a list of sh error files
ShFiles = [i for i in os.listdir() if '.sh.e' in i or '.sh.o' in i]

# loop of sh files
for filename in ShFiles:
    # get individual ID
    individual = filename[:filename.index('_')] 
    # extract datatype
    if 'RNA' in filename:
        datatype = 'RNAseq'
    elif 'GRCH37' in filename:
        datatype = 'WGS'
    # find the destination folder
    for tumor in MapTumor[datatype]:
        TumorFound = False
        if individual in MapTumor[datatype][tumor]:
            # tumor found, exit loop
            TumorFound = True            
            break
    
    # check that tumor has been found
    assert TumorFound == True, 'individual should be mapped to tumor'
    # move filename to appropriate folder   
    # check if folder already exists
    if datatype == 'RNAseq':
        destination = tumor + '/' + 'sh_files_rnaseq'
    elif datatype == 'WGS':
        destination = tumor + '/' + 'sh_files_wgs'
    try:
        os.listdir(destination)
    except:
        # create folder and move files
        os.mkdir(destination)
        os.system('mv ./' + filename + ' ' + destination)
    else:
        # folder already exists, move file to folder
        os.system('mv ./' + filename + ' ' + destination)
        
