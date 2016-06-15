# -*- coding: utf-8 -*-
"""
Created on Mon May  2 12:52:51 2016

@author: RJovelin
"""

# use this script to plot the CDF of allele frequencies for various functional categories

# place this script in folder with heteroplasmy summary files

# usage python3 PlotAlleleFrequencyFunction.py [options]
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
    
# create a dict {mutation: [list of allele frequencies]}
mutations = {}

# loop over filename in files
for filename in files:
    # compute allele frequencies for each mutation category
    # create a dict {functional_effect: [list of allele frequencies]}
    snpeffect = ComputeAlleleFrequency(filename)
    # pool all allelic frequencies togather
    for effect in snpeffect:
        # copy list of frequencies
        freq = list(snpeffect[effect])
        # add poitions to list value
        if effect in mutations:
            mutations[effect].extend(freq)
        else:
            mutations[effect] = freq


# sort frequency values
for effect in mutations:
    mutations[effect].sort()

# make parallel lists with functional effects and frequency values
categories = [i for i in mutations]
data = [mutations[i] for i in categories]

# replace category names, keeping the same orders between names and data
for i in range(len(categories)):
    if categories[i] == 'stopgain':
        categories[i] = 'PSC'
    elif categories[i] == 'stoploss':
        categories[i] = 'SCL'
    elif categories[i] == 'non-synonymous':
        categories[i] = 'NonSyn'
    elif categories[i] == 'synonymous':
        categories[i] = 'Syn'
    # keep the same name for the other categories (tRNA, DLoop, NonCoding, Ribosomal)

print(len(data), len(categories))
for i in categories:
    print(i)

# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

graph1 = ax.step(data[0], np.linspace(0, 1, len(data[0]), endpoint=False), linewidth = 1.2, color = '#8dd3c7', alpha = 0.7)
graph2 = ax.step(data[1], np.linspace(0, 1, len(data[1]), endpoint=False), linewidth = 1.2, color = '#ffffb3', alpha = 0.7)
graph3 = ax.step(data[2], np.linspace(0, 1, len(data[2]), endpoint=False), linewidth = 1.2, color = '#bebada', alpha = 0.7)
graph4 = ax.step(data[3], np.linspace(0, 1, len(data[3]), endpoint=False), linewidth = 1.2, color = '#fb8072', alpha = 0.7)
graph5 = ax.step(data[4], np.linspace(0, 1, len(data[4]), endpoint=False), linewidth = 1.2, color = '#80b1d3', alpha = 0.7)
graph6 = ax.step(data[5], np.linspace(0, 1, len(data[5]), endpoint=False), linewidth = 1.2, color = '#fdb462', alpha = 0.7)
graph7 = ax.step(data[6], np.linspace(0, 1, len(data[6]), endpoint=False), linewidth = 1.2, color = '#b3de69', alpha = 0.7)
graph8 = ax.step(data[7], np.linspace(0, 1, len(data[7]), endpoint=False), linewidth = 1.2, color = '#fccde5', alpha = 0.7)

print('plotted CDF')

# add label for the Y axis
ax.set_ylabel('Cumulative fraction of mutations', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Allele Frequency', size = 10, ha = 'center', fontname = 'Arial')

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
lns = graph1+graph2+graph3+graph4+graph5+graph6+graph7+graph8
# get labels
labs = [categories[i] for i in range(len(categories))]
# plot legend
ax.legend(lns, labs, loc=2, fontsize = 8, frameon = False)

# build outputfile namewith parameters
# extract the cancer name
if which_files == 'singlefile' and sample == 'tumor':
    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_' + sample)]
elif which_files == 'singlefile' and sample == 'specific':
    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_TumorSpecific')]

if which_files == 'singlefile':
    outputfile = 'CDFFreqMutations' + cancer + sample.capitalize() + '.pdf'
elif which_files == 'allfiles':
    outputfile = 'CDFFreqMutations' + sample.capitalize() + '.pdf'

fig.savefig(outputfile, bbox_inches = 'tight')