# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 15:40:08 2016

@author: RJovelin
"""

import os
import sys

print('done importing modules')

# use this script to prepare bash files and launch mitoseek on thse files
# Preconditin, place this script in /RQusagers/rjovelin/awadalla_group/TCGA/Mitoseek/ + (DestinationFolder) 
# use this script for unpaired bams and non-recalibrate WGS bams


# usage RunMitoSeekNonRecalibratedBams.py [options]
# - ID_file: file with participant ID and correspnding bam ID

# read bam info from file with participant ID and bam ID
IDFile = sys.argv[1]
print(IDFile)

# get the phred score from command
Phred = int(sys.argv[2])

print(IDFile)
print(Phred)


if Phred == 13:
    DestinationFolder = 'mitoseek_hp1_mbq13/'
elif Phred == 20:
    DestinationFolder = 'mitoseek_hp1_mbq20/'

# assign t = 4 by default, unless data is RNASeq
t = 4

# open file for reading
infile = open(IDFile, 'r')
# loop over file
for line in infile:
    print(line)
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get participant ID
        participant = line[0]
        print(participant)
        # get ref genome (not that MT genome ref is rCRS for Hg19, GRCH37 and GRCH37lie)
        genome = line[1].upper()
        print(genome)
        # get info on data type (tumor, normal, etc)
        datatype = line[3]
        # use only GRCH37 ref
        assert genome == 'GRCH37', 'genome is not GRCH37'
        myfile = line[2] + '.woPCRDUP.PP.UM.bam'
        
        # open file for writing command
        filepath = '/exec5/GROUP/awadalla/awadalla/awadalla_group/TCGA/Mitoseek/' + DestinationFolder
                
        outputfile = './' + participant + '_' + datatype + '_' + genome + '.sh' 
        newfile = open(outputfile, 'w')
        newfile.write('cd /RQusagers/rjovelin/awadalla_group/Programs/MitoSeek-1.3/' + '\n')
        newfile.write('module load bioperl/1.6.1\n')
        newfile.write('module load perl_extra\n')
        newfile.write('perl mitoSeek_debug.pl -i  /exec5/GROUP/awadalla/awadalla/awadalla_group/TCGA/Mitoseek/' + DestinationFolder +
        myfile + ' -t ' + str(t) + ' -sb 0 -sp 1 -mbq ' + str(Phred) + ' -hp 1 -ha 2 -r rCRS -R rCRS -o /exec5/GROUP/awadalla/awadalla/awadalla_group/TCGA/Mitoseek/' + DestinationFolder +
        participant + '_' + datatype + '_' + genome + ' -samtools /home/apps/Logiciels/samtools/0.1.19/bin/samtools' + '\n')
        newfile.close()            
        MyCommand = 'qsub -l walltime=2:00:00,nodes=1:ppn=4 ' + filepath + participant + '_' + datatype + '_' + genome + '.sh'
        os.system(MyCommand)
        print(MyCommand)
        
infile.close()
