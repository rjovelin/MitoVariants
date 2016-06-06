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
        to_delete.add(i)
# remve files and non-cancer folders
for i in to_delete:
    folders.remove(i)
print('N cancers', len(cancers))        










 BlackListedIndividuals.txt
 GenerateMitoseekSummaryFile.py
  rCRS_genes_MT.text.txt
