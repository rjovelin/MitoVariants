# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:24:55 2016

@author: RJovelin
"""

# use this script to generate the heteroplasmy summary files for RNA and WGS for all cancers
# Note that it doesn't generate the files for RNA Only


# place this script in mitoseek_hp1_mbq20

import os

# create a list of cancer names
cancers = os.listdir()

# remove all files and folders that are not cancer
to_remove = set()
for i in cancers:
    try:
        os.listdir(i)
    except:
        to_remove.add(i)
    if i == '__pycache__':
        to_remove.add(i)
# remve files and non-cancer folders
for i in to_remove:
    cancers.remove(i)
print('N cancers', len(cancers))        
print(cancers)

EssentialFiles = ['BlackListedIndividuals.txt', 'GenerateMitoseekSummaryFile.py', 'rCRS_genes_MT.text.txt']

# loop over cancer type
for i in range(len(cancers)):
    # copy essential files to TP_RNASeq and TP_WGS for all cancers
    for j in EssentialFiles:
        os.system('cp ' + j + ' ' + cancers[i] + '/TP_RNASeq/')
        os.system('cp ' + j + ' ' + cancers[i] + '/TP_WGS/')
print('done copying files')

# loop over cancer type
for i in range(len(cancers)):
    # move to TP_RNASeq in Cancer folder
    os.chdir(i + '/TP_RNASeq/')
    # generate summary file
    os.system('python3 GenerateMitoseekSummaryFile.py ' + cancers[i] + ' tumor RNA BlackListedIndividuals.txt 1000 1000')
    # move to TP_WGS in cancer folder
    os.system('../TP_WGS/')
    os.system('python3 GenerateMitoseekSummaryFile.py ' + cancers[i] + ' tumor WGS BlackListedIndividuals.txt 1000 1000')
    # move back to mitoseek_hp1_mbq20
    os.chdir('../../')

print('done generating summary files')
    