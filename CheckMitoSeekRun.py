# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:03:31 2016

@author: RJovelin
"""


# use this script to verify that MitoSeek ran correctly by reading the
# content of the bash error files

# usage: python3 CheckMitoSeekRun.py
# place this script in the same folder with thr bash error files

import os

files = [i for i in os.listdir() if '.sh.e' in i]
# go through files
for i in files:
    # open file for reading
    infile = open(i, 'r')
    assert '.bashrc' in infile.readline()
    assert '.bashrc' in infile.readline()
    line = infile.readline().rstrip()
    if line != '':
        print(i)
    infile.close()