# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:10:55 2016

@author: Richard
"""


# use this script to plot the proportions of SNPs with > 2 alleles for all tumors

# place this script in folder with heteroplasmy summary files

# usage python3 PlotMultiAlleles.py [options]
# - threshold: % reads to identifiy heteroplasmies
# - [frequency/counts]: plot frequency or counts of mutational effects
# - outputfile

import os
import sys
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
import numpy as np
from mito_mutations import *


# make a list of summary files 
files = [i for i in os.listdir() if 'RNAOnly' in i and '.txt' in i]

# get threshold from command
#threshold = int(sys.argv[1])
#frequency = sys.argv[2]
#outputfile = sys.argv[3]

outputfile = 'testfig.pdf'
frequency = 'counts'
threshold = 1

# create a dict {tumor: N multiallelic}
multialleles = {}

# loop over filename in files
for filename in files:
    # create a dict {mutation: {individuals: [list of positions]}}
    snps = IdentifyMultiAllelicPositions(filename, threshold)
    # count number of multiallele sites
    nums = 0
    for individual in snps:
        for site in snps[individual]:
            if len(snps[individual][site]) > 2:
                nums += 1
    # get tumor for filename
    tumor = filename[filename.index('_') + 1 : filename.index('_', filename.index('_') + 1)]
    multialleles[tumor] = nums    
       
# make a list of tumor names sorted by count
name_counts = []
for i in multialleles:
    name_counts.append([multialleles[i], i])
name_counts.sort()
tumors = [i[1] for i in name_counts]
counts = [i[0] for i in name_counts]
print(tumors)
print(counts)

# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add axe to fig
ax = fig.add_subplot(1, 1, 1)

# Set the bar width
bar_width = 0.4

# set positions of the left bar-boundaries
bar_left = [i for i in range(len(counts))]
# set positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+(bar_width/2) for i in bar_left]

# Create a bar plot, in position bar_left for counts
plt.bar(bar_left, counts, width=bar_width, color= 'red')

# set the x ticks with names
plt.xticks(tick_pos, tumors, size = 12)

#### need to edit the y ticks

if frequency == 'frequency':
    # count total number of mutations
    total = sum(counts)
    print(total)
    # set the y ticks
    counts = list(map(lambda x: x / total, counts))
    plt.yticks([i/100 for i in range(0, 125, 25)], [0, 0.25, 0.50, 0.75, 1])

# set axis labels
if frequency == 'frequency':
    plt.ylabel('Proportion of multiallelic mutations', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')
elif frequency == 'counts':
    plt.ylabel('Number of multiallelic mutations', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')

plt.xlabel('Tumor types', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# Set a buffer around the edge
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)


plt.margins()
  
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      

plt.title('Multiallelic sites in tumor genomes\n', size = 12)

plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black'
    )  

# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    