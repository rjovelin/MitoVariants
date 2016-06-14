# -*- coding: utf-8 -*-
"""
Created on Mon May  2 12:52:51 2016

@author: RJovelin
"""

# use this script to plot the CDF of allele frequencies for various functional categories

# place this script in folder with heteroplasmy summary files

# usage python3 PlotAlleleFrequencyFunction.py [options]
# - [singlefile/allfiles]: whether a single summary file or multiple summary files
# - [frequency/counts]: plot frequency or counts of mutational effects
# - [tumor/specific]: whether RNA variants are in tumor or are tumor specific (filtered based on normale)
# - [mutations/sites] : count the number of mutations or the number of variable sites
# - HeteroplasmySummaryFile: summary file of tumor if singlefile is used

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
# import custom modules
from mito_mutations import *

#use single of all summary files
which_files = sys.argv[1]
datatype = sys.argv[2] 
sample = sys.argv[3]
Sites = sys.argv[4]

if sample == 'tumor':
    # make a list of summary files 
    files = [i for i in os.listdir() if 'tumor_RNAOnly' in i and '.txt' in i]
elif sample == 'specific':
    # make a list of summary files 
    files = [i for i in os.listdir() if 'TumorSpecific' in i and '.txt' in i]
    
if which_files == 'singlefile':
    # get tumor from command
    HeteroplasmySummaryFile = sys.argv[5]
    files = [filename]



#################### CONTINUE HERE





# create a dict {mutation: [list of positions]}
mutations = {}

# loop over filename in files
for filename in files:
    # create a dict {mutation: [list of positions]}
    snpeffect = ComputeMAFFunction(filename)
    # populate mutation dict
    for effect in snpeffect:
        # copy list of frequencies
        freq = list(snpeffect[effect])
        # add poitions to list value
        if effect in mutations:
            mutations[effect].extend(freq)
        else:
            mutations[effect] = freq

categories = [i for i in mutations if i != 'NA']
print(categories)
snps = []
for i in mutations:
    if i != 'NA':
        snps.append(mutations[i])
plt.hist(snps, label=categories)        
        
#plt.hist(snps, label = categories)

plt.show()

var = {}
for i in mutations:
    if i != 'NA':
        for j in mutations[i]:
            if j < 0.1:
                if i in var:
                    var[i].append(j)
                else:
                    var[i] = [j]
muts = []
for i in var:
    muts.append([len(var[i]), i])
muts.sort()
print(muts)


# add a title
# add a label
# order freqs in decreasing order of the first freq bin
# change y and x ticks
# opetion to have frequencies instead of counts


            
## make a list of mutation names sorted by count
#name_counts = []
#for i in mutations:
#    if i in ['DLoop', 'non-synonymous', 'stoploss', 'stopgain', 'synonymous', 'NonCoding', 'tRNA']:
#        if i == 'non-synonymous':
#            name_counts.append([len(mutations[i]), 'NonSyn'])
#        elif i == 'synonymous':
#            name_counts.append([len(mutations[i]), 'Syn'])
#        elif i == 'stoploss':
#            name_counts.append([len(mutations[i]), 'SCL'])
#        elif i == 'stopgain':
#            name_counts.append([len(mutations[i]), 'PSC'])
#        else:
#            name_counts.append([len(mutations[i]), i])
#name_counts.sort()
#categories = []
#for i in name_counts:
#    categories.append(i[1])
#    
#counts = [i[0] for i in name_counts]
#print(categories)
#print(counts)
#
## create figure
#fig = plt.figure(1, figsize = (4.3,2.56))
#
## add axe to fig
#ax = fig.add_subplot(1, 1, 1)
#
## Set the bar width
#bar_width = 0.4
#
## set positions of the left bar-boundaries
#bar_left = [i for i in range(len(counts))]
## set positions of the x-axis ticks (center of the bars as bar labels)
#tick_pos = [i+(bar_width/2) for i in bar_left]
#
## Create a bar plot, in position bar_left for counts
#plt.bar(bar_left, counts, width=bar_width, color= 'red')
#
## set the x ticks with names
#plt.xticks(tick_pos, categories, size = 12)
#
##### need to edit the y ticks
#
#if frequency == 'frequency':
#    # count total number of mutations
#    total = sum(counts)
#    print(total)
#    # set the y ticks
#    counts = list(map(lambda x: x / total, counts))
#    plt.yticks([i/100 for i in range(0, 125, 25)], [0, 0.25, 0.50, 0.75, 1])
#
## set axis labels
#if frequency == 'frequency':
#    plt.ylabel('Proportion of mutations', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')
#elif frequency == 'counts':
#    plt.ylabel('Number of mutations', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')
#
#plt.xlabel('Mutational effect', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
#
## Set a buffer around the edge
#plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
#
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
## hide these grids behind plot objects
#ax.set_axisbelow(True)
#
#
#plt.margins()
#  
## do not show lines around figure  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(False)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(False)      
#
#if which_files == 'singlefile':
#    plt.title(tumor, size = 12)
#elif which_files == 'allfiles':
#    plt.title('All tumor types', size = 12)
#
#plt.tick_params(
#    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='on', # labels along the bottom edge are off 
#    colors = 'black'
#    )  

# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    
    
    
    
    
    
    
    
    
##########################
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 07:21:34 2016
@author: Richard
"""



# use this script to plot the cumulative distribution function of the MAF distribution for miRNAs and other sites

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import sys
import json
# import custom modules
from manipulate_sequences import *
from miRNA_target import *
from genomic_coordinates import *
from sites_with_coverage import *
from divergence import *
from premature_stops import *
from randomize_SNPs import *
from get_coding_sequences import *


# make a list with MAF of replacement SNPs, excluding sites with sample size < 10
MAF_REP = MAF_SNP('../Genome_Files/CDS_SNP_DIVERG.txt', '../Genome_Files/unique_transcripts.txt', 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'REP', 10)
print('REP', len(MAF_REP))
print('MAF for replacement sites done')

# make a list with MAF of synonymous sites, excluding sites with sample size < 10
MAF_SYN = MAF_SNP('../Genome_Files/CDS_SNP_DIVERG.txt', '../Genome_Files/unique_transcripts.txt', 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'SYN', 10)
print('SNP', len(MAF_SYN))
print('MAF for synonymous sites done')

# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts in genome')

# get the coordinates of the miRNA loci
mirnas_coord = get_mirna_loci('CRM_miRNAsCoordinatesFinal.txt')
print('got miRNA coordinates')
# get the allele counts for miRNA sites
mirna_sites = get_feature_sites(chromo_sites, mirnas_coord)
print('got allele counts for miRNAs')
# compute MAF for mirna sites (sites with sample size < 10 are already excluded)
MAF_mirna = MAF_non_coding(mirna_sites)
print('miRNAs', len(MAF_mirna))
print('MAF for miRNA sites done')

# get target coordinates from json file
infile = open('AllmiRNATargetsCoords.json')
target_coord = json.load(infile)
infile.close()
print('target coords', len(target_coord))
print('got miRNA target site coordinates')
# get the allele counts for all targets
target_sites = get_feature_sites(chromo_sites, target_coord)
print('got allele counts for target sites')
# compute MAF for all targets
MAF_targets = MAF_non_coding(target_sites)
print('MAF for miRNA targets done')

# sort MAF values
MAF_REP.sort()
MAF_SYN.sort()
MAF_mirna.sort()
MAF_targets.sort()
print('values are sorted')

# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot MAF synonymous sites
graph1 = ax.step(MAF_SYN, np.linspace(0, 1, len(MAF_SYN), endpoint=False), linewidth = 1.2, color = '#33a02c', alpha = 0.7)
# plot MAF replacement sites
graph2 = ax.step(MAF_REP, np.linspace(0, 1, len(MAF_REP), endpoint=False), linewidth = 1.2, color = '#b2df8a', alpha = 0.7)
# plot MAF miRNAs
graph3 = ax.step(MAF_mirna, np.linspace(0, 1, len(MAF_mirna), endpoint=False), linewidth = 1.2, color = '#1f78b4', alpha = 0.7)
# plot MAF targets
graph4 = ax.step(MAF_targets, np.linspace(0, 1, len(MAF_targets), endpoint=False), linewidth = 1.2, color = '#a6cee3', alpha = 0.7)
print('plotted CDF')

# add label for the Y axis
ax.set_ylabel('Proportion of SNPs', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Minor Allele Frequency', size = 10, ha = 'center', fontname = 'Arial')

# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# do not show ticks on 1st graph
ax.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

# do not show ticks
ax.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# add lines
lns = graph1+graph2+graph3+graph4
# get labels
labs = ['Synonymous', 'Replacement', 'miRNAs', 'targets']
# plot legend
ax.legend(lns, labs, loc=2, fontsize = 8, frameon = False)

fig.savefig('CDFMAFmiRNAs.pdf', bbox_inches = 'tight')