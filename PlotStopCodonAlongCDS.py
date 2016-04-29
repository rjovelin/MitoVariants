# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:28:04 2016

@author: RJovelin
"""

import os
import sys
from mito_mutations import *
import matplotlib.pyplot as plt


# create an histogramm with the number of mutations creating a stop codon along the CDS

# place this script in the same folder with heteroplasmy summary files

# usage python3 PlotStopCodonAlongCDS.py [options]
# - [RNA/Tumor]: RNA variants only or tumor specific variants
# outputfile: save figure to outputfile

datatype = sys.argv[1]
outputfile = sys.argv[2]

# create a list of heteroplasmy summary files
if datatype == 'RNA':
    # record variants found only in RNA
    files = [i for i in os.listdir() if 'RNAOnly' in i and '.txt' in i]
elif datatype == 'Tumor':
    files = [i for i in os.listdir() if 'Specific' in i and '.txt' in i]

# create a list to count all stop codons mutations
PTC = []
# loop over file name, record stop codon mutations
for filename in files:
    stops = StopCodonAlongCDS(filename, 'rCRS_genes_MT.text.txt', Sites = 'mutations', FirstOnly = False)
    PTC.extend(stops)


# create figure
fig = plt.figure(1, figsize = (4.3,2.56))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# create histogram
ax.hist(PTC, range(0, 110, 10), color = [0.392, 0.647, 0.769], edgecolor = 'w')

# add title
ax.set_title('Stop codon mutations along coding sequences\n', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# set y axis label
ax.set_ylabel('Number of PSC mutations', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# add labels to x-ticks, rotate and align right
ax.set_xticklabels(range(0, 110, 10), rotation = 0, ha = 'center', size = 10, fontname = 'Helvetica', family = 'sans-serif')

plt.yticks(fontsize = 10)
plt.xticks(range(0, 110, 10))

# set x axis label
ax.set_xlabel('Relative CDS length (%)', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)
# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.margins()
  
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
  
# do not show ticks
plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are on
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  
  
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')


