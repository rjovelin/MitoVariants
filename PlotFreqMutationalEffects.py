# -*- coding: utf-8 -*-
"""
Created on Sun May  1 08:02:51 2016

@author: Richard
"""

# use this script to plot the proportions of SNPs in various functional categories

# place this script in folder with heteroplasmy summary files

# usage python3 PlotFreqFunctionalEffect.py [options]
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

# create a dict {mutation: [list of positions]}
mutations = {}

# loop over filename in files
for filename in files:
    # create a dict {mutation: [list of positions]}
    snpeffect = MutationalEffect(filename, Sites)
    # populate mutation dict
    for effect in snpeffect:
        # copy list of positions
        positions = list(snpeffect[effect])
        # add poitions to list value
        if effect in mutations:
            mutations[effect].extend(positions)
        else:
            mutations[effect] = positions
            
for i in mutations:
    print(i, len(mutations[i]))            
            
          
# make a list of mutation names sorted by count
name_counts = []
for i in mutations:
    if i in ['DLoop', 'non-synonymous', 'stoploss', 'stopgain', 'synonymous', 'NonCoding', 'tRNA', 'Ribosomal']:
        if i == 'non-synonymous':
            name_counts.append([len(mutations[i]), 'NonSyn'])
        elif i == 'synonymous':
            name_counts.append([len(mutations[i]), 'Syn'])
        elif i == 'stoploss':
            name_counts.append([len(mutations[i]), 'SCL'])
        elif i == 'stopgain':
            name_counts.append([len(mutations[i]), 'PSC'])
        else:
            name_counts.append([len(mutations[i]), i])
name_counts.sort()
categories = []
for i in name_counts:
    categories.append(i[1])
    
counts = [i[0] for i in name_counts]
print(categories)
print(counts)

colorscheme = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5']

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
plt.bar(bar_left, counts, width=bar_width, color= colorscheme)

# set the x ticks with names
plt.xticks(tick_pos, categories, ha = 'right', rotation = 20, size = 12)


if datatype == 'frequency':
    # count total number of mutations
    total = sum(counts)
    print(total)
    # set the y ticks
    counts = list(map(lambda x: x / total, counts))
    #plt.yticks([i/100 for i in range(0, 125, 25)], [0, 0.25, 0.50, 0.75, 1])

# set axis labels
if datatype == 'frequency':
    plt.ylabel('Proportion of mutations', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')
elif datatype == 'counts':
    plt.ylabel('Number of mutations', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')

plt.xlabel('Mutational effect', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

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

if which_files == 'singlefile':
    if sample == 'tumor':
        plt.title('RNA variants in {0}'.format(tumor), size = 12)
    elif sample == 'specific':
        plt.title('Tumor-specific RNA variants in {0}'.format(tumor), size = 12)
elif which_files == 'allfiles':
    if sample == 'tumor':
        plt.title('RNA variants across all tumor types', size = 12)
    elif sample == 'specific':
        plt.title('Tumor-specific RNA variants in LIRI and RECA', size = 12)

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

if which_files == 'singlefile':
    # extract the cancer name
    if sample == 'tumor':
        cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_' + sample)]
    elif sample == 'specific':
        cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_TumorSpecific')]
           
# build outputfile with parameters  
if which_files == 'singlefile':
    outputfile = 'MutationalEffect' + cancer + sample.capitalize() + datatype.capitalize() + '.pdf' 
elif which_files == 'allfiles':
    outputfile = 'MutationalEffect' + sample.capitalize() + datatype.capitalize() + Sites.capitalize() + '.pdf' 
    
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')    