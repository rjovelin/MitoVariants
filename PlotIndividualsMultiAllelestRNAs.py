# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 11:49:41 2016

@author: RJovelin
"""


# create a figure with number of individuals with multiallelic sites for each positions of tRNAs

# place this script in folder with heteroplasmy summary files

# usage python3 PlotIndividualsMultiAllelestRNAs.py [options]
# - threshold: % reads to identifiy heteroplasmies
# - [specific/tumor]: consider tumor specific RNA variants (variable positions in normal are filtered)
#                     or tumor RNA variants (include both tissue specific and tumor specific variants)
# [singlefile, allfiles]: if a single summary file is used or if data is pooled across all summary files

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


# get threshold from command
threshold = float((sys.argv[1]))
sample = sys.argv[2]
which_files = sys.argv[3]

# make a list of summary files 
if sample == 'tumor' and which_files == 'allfiles':
    files = [i for i in os.listdir() if 'RNAOnly' in i and '.txt' in i]
elif sample == 'specific' and which_files == 'allfiles':
    files = [i for i in os.listdir() if 'TumorSpecific' in i and '.txt' in i]
elif which_files == 'singlefile':
    # get tumor from command
    HeteroplasmySummaryFile = sys.argv[4]
    files = [HeteroplasmySummaryFile]
    # get cancer name from summary file
    if sample == 'specific':
        cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_TumorSpecific')]
    else:
        cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_' + sample)]

# get the positions of each mitochiondrial gene
mito_genes = MitoAnnotation('rCRS_genes_MT.text.txt')
# make a set of tRNA positions
trna_indices = []
for gene in mito_genes:
    if gene.startswith('TRN'):
        for i in mito_genes[gene]:
            trna_indices.append(i)

# get the gene coordinates [start, end, orientation]
mito_coords = MitoCoordinates('rCRS_genes_MT.text.txt')

# create a dict {positions: N individuals with multiallelic variant}
multialleles = {}

# loop over filename in files
for filename in files:
    # create a dict {individuals: {positions: [alleles]}}
    snps = IdentifyMultiAllelicPositions(filename, threshold)
    # check if variant site is multiallelic
    for individual in snps:
        for site in snps[individual]:
            # check if site in tRNA
            if site in trna_indices:
                # check if variant is multiallelic
                if len(snps[individual][site]) > 2:
                    # find the tRNA corresponding to current site
                    for gene in mito_genes:
                        if site in mito_genes[gene]:
                            # stop looking, found tRNA gene
                            break
                    # convert genomic position to tRNA position
                    assert site in range(mito_coords[gene][0], mito_coords[gene][1]), 'position should be in {0}'.format(gene)
                    position = GenomicPositionToGenePosition(site, mito_coords[gene][0], mito_coords[gene][1], mito_coords[gene][2])                    
                    # populate dict
                    if position in multialleles:
                        multialleles[position] += 1
                    else:
                        multialleles[position] = 1
        
# get a list of indices corresponding to the longest tRNAs
# (ie. all tRNAs are aligned at first nt to plot N at each position)
longest = 0
for gene in mito_coords:
    if gene.startswith('TRN'):
        # get the length of the tRNA
        L = mito_coords[gene][1] - mito_coords[gene][0] + 1
        # update longest
        if L > longest:
            longest = L

# add poisitions in dict if not present
for i in range(longest):
    if i not in multialleles:
        multialleles[i] = 0

# create a list of positions
positions = [i for i in multialleles]
positions.sort()
print('tRNA ', len(positions), min(positions), max(positions))
# create a list of counts
counts = [multialleles[i] for i in positions]

# create parallel lists with counts and positions, without 0s
Indices = [i for i in range(len(counts)) if counts[i] != 0]
Counts = [counts[i] for i in Indices]
# create a list of positions 1-index based
Positions = [positions[i]+1 for i in Indices]
print(Counts)
print(Positions)

# create figure
fig = plt.figure(1, figsize = (3.5, 2))
# add axe to fig
ax = fig.add_subplot(1, 1, 1)

# Set the bar width
bar_width = 0.2

# set up color
colorscheme = ['red', 'black', 'black']

# Create a bar plot, in position bar_left for counts
for i in range(len(Positions)):
    plt.bar(Positions[i], Counts[i], width = bar_width, color= colorscheme[i])

# set positions of the x-axis ticks (center of the bars as bar labels)
#tick_pos = [i+(bar_width/2) for i in bar_left]
# set the x ticks with names
#plt.xticks(tick_pos, tumors, rotation = 20, ha = 'right', size = 12)

# limit the x axis
plt.xlim([0, max(positions)+1])

# add axes labels
plt.ylabel('Number of individuals', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'black')
plt.xlabel('Positions in tRNAs', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# add margins
plt.margins(0.05)
  
# do not show lines around figure  
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


plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    direction = 'out'
    )  


if which_files == 'singlefile':
    outputfile = 'MultiAlleleCountstRNAPositions' + cancer + sample.capitalize() + 'Het' + str(threshold) + '.pdf'
elif which_files == 'allfiles':
    outputfile = 'MultiAlleleCountstRNAPositions' + sample.capitalize() + 'Het' + str(threshold) + '.pdf'

# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    
