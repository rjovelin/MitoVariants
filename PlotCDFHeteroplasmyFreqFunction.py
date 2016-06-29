# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 15:23:39 2016

@author: RJovelin
"""

# use this script to plot the CDF of heteroplasmie frequencies (ie. within individual, a same mutant may be recorded multiple times)
# for various functional categories


# place this script in folder with heteroplasmy summary files

# usage python3 PlotCDFHeteroplasmyfreqFunction.py [options]
# - [singlefile/allfiles]: whether a single summary file or multiple summary files
# - [tumor/specific]: whether RNA variants are in tumor or are tumor specific (filtered based on normale)
# - threshold (in%) to detect heteroplasmy
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
threshold = float(sys.argv[3])

if which_files == 'allfiles' and sample == 'tumor':
    # make a list of summary files 
    files = [i for i in os.listdir() if 'tumor_RNAOnly' in i and '.txt' in i]
elif which_files == 'allfiles' and sample == 'specific':
    # make a list of summary files 
    files = [i for i in os.listdir() if 'TumorSpecific' in i and '.txt' in i]
elif which_files == 'singlefile':
    # get tumor from command
    files = [sys.argv[4]]
    if sample == 'tumor':
        assert 'tumor_RNAOnly' in files[0], 'summary file should be for RNA variants'
    elif sample == 'specific':
        assert 'TumorSpecific' in files[0], 'summary file sould be for tumor-specific RNA variants'
    
 
# create a dict {mutation: [list of allele frequencies]}
mutations = {}

# loop over filename in files
for filename in files:
    # compute allele frequencies for each mutation category
    # create a dict {functional_effect: [list of mutant allele frequencies]}
    snpeffect = ComputeHeteroplasmyFrequencyMutant(filename, threshold)
    # pool all allelic frequencies togather
    for effect in snpeffect:
        # add frequencies to list value
        if effect in mutations:
            mutations[effect].extend(list(snpeffect[effect]))
        else:
            mutations[effect] = list(snpeffect[effect])

# sort frequency values
for effect in mutations:
    mutations[effect].sort()

# make list with functional category and data
Data = [[key, val] for key, val in mutations.items()]
# replace some functional categories with a shorter name
for i in range(len(Data)):
    if Data[i][0] == 'stopgain':
        Data[i][0] = 'PSC'
    elif Data[i][0] == 'stoploss':
        Data[i][0] = 'SCL'
    elif Data[i][0] == 'non-synonymous':
        Data[i][0] = 'NonSyn'
    elif Data[i][0] == 'synonymous':
        Data[i][0] = 'Syn'
    # keep the same name for the other categories (tRNA, DLoop, NonCoding, Ribosomal)

# sort according to names
Data.sort()

# create parallel lists with functional names and data
categories = [Data[i][0] for i in range(len(Data))]
data = [Data[i][1] for i in range(len(Data))]    
print(categories)


# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# create a list to store the variables for building the legend
Graphs = []

# make a list of colors
colorscheme = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5']

# loop over categories
for i in range(len(categories)):
    print(categories[i])
    graph = ax.step(data[i], np.linspace(0, 1, len(data[i]), endpoint=True), linewidth = 1.2, color = colorscheme[i], alpha = 0.7)
    Graphs.append(graph)
print('plotted CDF')

# add label for the Y axis
ax.set_ylabel('Cumulative fraction of mutations', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Mutant allele Frequency', size = 10, ha = 'center', fontname = 'Arial')

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
lns = Graphs[0]
for i in range(1, len(Graphs)):
    lns = lns + Graphs[i]
# get labels
labs = [categories[i] for i in range(len(categories))]
print(labs)
# plot legend
ax.legend(lns, labs, loc=4, fontsize = 6, frameon = False)

## build outputfile namewith parameters
## extract the cancer name
#if which_files == 'singlefile' and sample == 'tumor':
#    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_' + sample)]
#elif which_files == 'singlefile' and sample == 'specific':
#    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_TumorSpecific')]
#
#if which_files == 'singlefile':
#    outputfile = 'CDFFreqMutations' + cancer + sample.capitalize() + '.pdf'
#elif which_files == 'allfiles':
#    outputfile = 'CDFFreqMutations' + sample.capitalize() + '.pdf'

fig.savefig('testfile.pdf', bbox_inches = 'tight')