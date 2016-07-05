# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 11:54:09 2016

@author: RJovelin
"""

# use this script to plot correlation between genotype at MRPP3 and RNA heteroplasmy frequency

# usage PlotCorrelMRPP3GenotypeAltFreq [options]
# -[trna/P9/allgenes]: consider only RNA modifications in tRNA, modifs in tRNA P9 or all genes


# place this script in folder HeteroPlasmyFilesRNAOnly with summary files

# import matplotlib and change api to use on server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import built in modules
import sys
import os
import numpy as np
from scipy import stats
# import custom modules
from mito_mutations import *



# create a dict {individual: genotype}
Genotypes = {}

# make a list of tumor types with VCF files
VCFTumor = ['COAD', 'OV', 'RECA-EU', 'UCEC', 'CESC', 'LGG', 'LIRI', 'SARC', 'STAD']

# loop over tumor folders, unzip VCF files
for folder in VCFTumor:
    # make a list of gzip files    
    if folder == 'STAD':
        files = [i for i in os.listdir(folder + '/filteredVCFS') if i[-3:] == '.gz']
    else:
        files = [i for i in os.listdir(folder) if i[-3:] == '.gz']
    if len(files) != 0:
        # unzip files
        if folder == 'STAD':
            for filename in files:
                os.system('gunzip ' + folder + '/filteredVCFS/' + filename)
        else:
            for filename in files:
                os.system('gunzip ' + folder + '/' + filename)
print('done unzipping VCF files')                

# make a list of WGS normal files with corresponding individual ID and 
MatchFiles = [i for i in os.listdir() if 'WGS_NTOnly' in i]
assert len(MatchFiles) == 9, 'there should be matching files for 9 tumors'

# create a dict {BAM ID : individual ID} 
MatchingIDs = {}
for filename in MatchFiles:
    infile = open(filename)
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            # get individual ID
            Individual = line[0]
            # get BAM ID
            BamID = line[2]
            BamID = BamID[BamID.index('wgsmito/') + len('wgsmito/'): BamID.index('.mito.')]
            # populate dict
            MatchingIDs[BamID] = Individual
    infile.close()
print('matched Bams with Individual IDs')



 


           
            


# find the ID - VCF matchs

# get the genotypes of MRPP3 at rs11156878 

# get the heteroplasmy frequency for trna, p9 or all genes 





/RQusagers/rjovelin/awadalla_group/TCGA/VCFs/Germlines_vcf/perCancerType


/RQusagers/rjovelin/awadalla_group/TCGA/Mitoseek/mitoseek_hp1_mbq20

/RQusagers/rjovelin/awadalla_group/TCGA/Mitoseek/mitoseek_hp1_mbq20/HeteroPlasmyFilesRNAOnly




# create a correlation plot between alternative allele frequency as a function of genotypes at given position

# usage plot_correlation_altfreq_genotype.py [parameters]
# inputfile : file with genotypes and allele frequencies
# outputfile : figure file
# SNP ID: ID of the SNP with genotype






# get inputfile from command
genotype_file = sys.argv[1]
# get outputfile
figure_file = sys.argv[2]
# get the ID of the SNP
SNP_ID = sys.argv[3]

# create a dict with {genotypes: list of allelic frequencies}
genotypes = {}
# open file for reading
infile = open(genotype_file, 'r')
#loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get genotype
        genotype = line[1]
        # make a list of allelic frequencies
        freq = line[2].split(';')
        altfreq = list(map(lambda x: float(x), freq))
        # check if genotype in dict
        if genotype in genotypes:
            # add all frequency values
            genotypes[genotype].extend(altfreq)
        else:
            # add genotype : list of frequencies pairs to dict
            genotypes[genotype] = altfreq
# close file 
infile.close()


# make a list of genotypes
SNP = [i for i in genotypes]
# sort list
SNP.sort()


# create a figure
fig = plt.figure(1, figsize = (5, 3))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)

#for i in range(len(SNP)):
#    for j in genotypes[SNP[i]]:
#        ax.scatter(i+1, j, edgecolor = 'black', facecolor = 'black', lw = 1, s = 10, alpha = 0.5)



xA = np.random.normal(1, 0.1, len(genotypes[SNP[0]])) 
         

xB = np.random.normal(2, 0.1, len(genotypes[SNP[1]]))

xC = np.random.normal(3, 0.1, len(genotypes[SNP[2]]))


plt.scatter(xA, genotypes[SNP[0]], edgecolor = 'black', facecolor = 'black', lw = 1, s = 10, alpha = 0.5)
plt.scatter(xB, genotypes[SNP[1]], edgecolor = 'grey', facecolor = 'grey', lw = 1, s = 10, alpha = 0.5) 
plt.scatter(xC, genotypes[SNP[2]], edgecolor = 'lightgrey', facecolor = 'lightgrey', lw = 1, s = 10, alpha = 0.5)


ax.set_ylim([0, 1])
ax.set_xlim([0,4])
            
# set y axis label
ax.set_ylabel('Alternative Allele Frequency', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels([SNP[0], SNP[1], SNP[2]], ha = 'center', size = 10, fontname = 'Helvetica', family = 'sans-serif')

# set x axis label
ax.set_xlabel(SNP_ID, size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


# set positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+1 for i in range(3)]

plt.xticks(tick_pos, SNP)














# save figure
fig.savefig(figure_file, bbox_inches = 'tight')




