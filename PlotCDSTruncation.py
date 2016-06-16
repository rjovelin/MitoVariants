# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 12:28:37 2016

@author: RJovelin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:28:04 2016

@author: RJovelin
"""


# create an histogramm with the number of mutations creating a stop codon along the CDS

# place this script in the same folder with heteroplasmy summary files

# usage python3 PlotStopCodonAlongCDS.py [options]
# - [singlefile/allfiles]: whether a single summary file or multiple summary files
# - [tumor/specific]: whether RNA variants are in tumor or are tumor specific (filtered based on normale)
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
sample = sys.argv[2]

if which_files == 'allfiles' and sample == 'tumor':
    # make a list of summary files 
    files = [i for i in os.listdir() if 'tumor_RNAOnly' in i and '.txt' in i]
elif which_files == 'allfiles' and sample == 'specific':
    # make a list of summary files 
    files = [i for i in os.listdir() if 'TumorSpecific' in i and '.txt' in i]
elif which_files == 'singlefile':
    # get tumor from command
    files = [sys.argv[3]]
print(files)    


# create a list to count all stop codons mutations
PTC = []
# loop over file name, record stop codon mutations
for filename in files:
    stops = StopCodonAlongCDS(filename, 'rCRS_genes_MT.text.txt', Sites = 'mutations', FirstOnly = False)
    PTC.extend(stops)


# create figure
fig = plt.figure(1, figsize = (4, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2

# create histogram
ax.hist(PTC, range(0, 110, 10), color = '#9ecae1', edgecolor = '#9ecae1')

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


# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

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

# add margins around x axis 
plt.margins(0.05)


# build output file name from parameters

# extract the cancer name
if which_files == 'singlefile' and sample == 'tumor':
    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_' + sample)]
elif which_files == 'singlefile' and sample == 'specific':
    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_TumorSpecific')]

if which_files == 'singlefile':
    outputfile = 'StopCodonsDistribution' + cancer + sample.capitalize() + '.pdf'
elif which_files == 'allfiles':
    outputfile = 'StopCodonsDistribution' + sample.capitalize() + '.pdf'

# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
