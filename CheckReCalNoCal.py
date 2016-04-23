# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 01:52:46 2016

@author: Richard
"""

import os
import sys

# use this script to prepare bash files and launch mitoseek on thse files
# to compare heteroplasmies in same individual with calibrated and no calibrated bams

# Precondition, place this script in /RQusagers/rjovelin/awadalla_group/TCGA/Mitoseek/ + (DestinationFolder) 

# usage CheckReCalNoCal.py [options]
# - [nocal/recal]: use recalibrated or no calibrated bams

# get option from command
calibration = sys.argv[1]

# make a list of input bams
if calibration == 'recal':
    bams = [i for i in os.listdir() if i[-4:] == '.bam' and 'recalibrate' in i]
elif calibration == 'nocal':
    bams = [i for i in os.listdir() if i[-4:] == '.bam' and 'recalibrate' not in i]


# loop over bam files
for filename in bams:
    # get participant ID
    participant = filename[:filename.index('.')]
    print(participant)
    
    # open file for writing command
    filepath = '/exec5/GROUP/awadalla/awadalla/awadalla_group/TCGA/Mitoseek/norecalVSrecal/'
    outputfile = './' + participant + '_' + calibration + '_TP_GRCH37' + '.sh' 
    newfile = open(outputfile, 'w')
    newfile.write('cd /RQusagers/rjovelin/awadalla_group/Programs/MitoSeek-1.3/' + '\n')
    newfile.write('module load bioperl/1.6.1\n')
    newfile.write('module load perl_extra\n')
    newfile.write('perl mitoSeek_debug.pl -i  /exec5/GROUP/awadalla/awadalla/awadalla_group/TCGA/Mitoseek/norecalVSrecal/' + filename + 
    ' -t 4 -sb 0 -sp 1 -mbq 20 -hp 1 -ha 2 -r rCRS -R rCRS -o /exec5/GROUP/awadalla/awadalla/awadalla_group/TCGA/Mitoseek/norecalVSrecal/' +
    participant + '_TP_GRCH37 -samtools /home/apps/Logiciels/samtools/0.1.19/bin/samtools' + '\n')
    newfile.close()            
    MyCommand = 'qsub -l walltime=2:00:00,nodes=1:ppn=4 ' + filepath + participant + '_' + calibration + '_TP_GRCH37' + '.sh'
    os.system(MyCommand)
    print(MyCommand)










