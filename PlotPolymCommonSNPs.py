# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:30:46 2016

@author: RJovelin
"""



# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 15:18:49 2016

@author: RJovelin
"""

import matplotlib.pyplot as plt
import sys
from mito_mutations import *


# create a figure with polymorphism content at variable positions shared among tumor types

# place script in same folder with heteroplasmy summary files

# usage python3 PlotPolymCommonSNPs.py [parameters]
# - threshold: the value above which PIC is considered (eg. threshold = 0 means 0 values won't be shown)
# - [normal/cancer] whether the sample is germline or tumor 
# - [full/limited/black] : color panel: each gene has a different color (full) 
#                        or proteins are in grey and tRNAs in red (limited)
#                        or only tRNAs are labeled in red (black)
# - outputfile: save fig to outputfile


## get thresold from command
#threshold = float(sys.argv[1])
## get normal or cancer from command
#sample = sys.argv[2]
## get the color panel to use
#colors = sys.argv[3]
## get outputfile from command
#outputfile = sys.argv[6]


threshold = 0
# get normal or cancer from command
sample = 'cancer'
# get the color panel to use
colors = 'limited'
# get outputfile from command
outputfile = 'testfig.pdf'


# make a list of hetero summary files
files = [i for i in os.listdir() if 'RNAOnly' in i and '.txt' in i]

# create a set of common SNPs among tumor types
commonsnps = CommonVariablePositions(files)

# get the positions of each mitochiondrial genes and regions
mito_genes = MitoAnnotation('rCRS_genes_MT.text.txt')

# create a dict to with common snps : list of PIC values pairs
Polym = {}
for i in commonsnps:
    Polym[i] = []

# loop over summary file
for filename in files:
    # compute PIC for each individual at each position
    polymorphism = GenomePositionsPic(filename)
    # add PIV values to common snps
    for i in polymorphism:
        if i in Polym:
            Polym[i].extend(polymorphism[i])

# create a list of positions
positions = [i for i in Polym]

# check which color panel to use
if colors == 'full':
    # create a dict with gene : color
    gene_colors = {'ATP6': 'yellow', 'ATP8': 'blue', 'COX1': 'lightgrey', 'COX2': 'blue',
               'COX3': 'green', 'CYTB': 'blue', 'ND1': 'yellow', 'ND2': 'orange',
               'ND3': 'magenta', 'ND4': 'green', 'ND4L': 'magenta', 'ND5': 'pink',
               'ND6': 'grey', 'RNR1': 'green', 'RNR2': 'lightgrey', 'TRNA': 'red',
               'TRNC': 'red', 'TRND':'red', 'TRNE': 'red', 'TRNF': 'red', 'TRNG': 'red',
               'TRNH': 'red', 'TRNI': 'red', 'TRNK': 'red', 'TRNL1': 'red', 'TRNL2': 'red',
               'TRNM': 'red', 'TRNN': 'red', 'TRNP': 'red', 'TRNQ': 'red', 'TRNR': 'red',
               'TRNS1': 'red', 'TRNS2':'red', 'TRNT': 'red', 'TRNV': 'red', 'TRNW': 'red',
               'TRNY': 'red'}
elif colors == 'black':
    gene_colors = {'ATP6': 'black', 'ATP8': 'black', 'COX1': 'black', 'COX2': 'black',
               'COX3': 'black', 'CYTB': 'black', 'ND1': 'black', 'ND2': 'black',
               'ND3': 'black', 'ND4': 'black', 'ND4L': 'black', 'ND5': 'black',
               'ND6': 'black', 'RNR1': 'black', 'RNR2': 'black', 'TRNA': 'red',
               'TRNC': 'red', 'TRND':'red', 'TRNE': 'red', 'TRNF': 'red', 'TRNG': 'red',
               'TRNH': 'red', 'TRNI': 'red', 'TRNK': 'red', 'TRNL1': 'red', 'TRNL2': 'red',
               'TRNM': 'red', 'TRNN': 'red', 'TRNP': 'red', 'TRNQ': 'red', 'TRNR': 'red',
               'TRNS1': 'red', 'TRNS2':'red', 'TRNT': 'red', 'TRNV': 'red', 'TRNW': 'red',
               'TRNY': 'red'}

elif colors == 'limited':
    gene_colors = {'ATP6': 'grey', 'ATP8': 'grey', 'COX1': 'grey', 'COX2': 'grey',
               'COX3': 'grey', 'CYTB': 'grey', 'ND1': 'grey', 'ND2': 'grey',
               'ND3': 'grey', 'ND4': 'grey', 'ND4L': 'grey', 'ND5': 'grey',
               'ND6': 'grey', 'RNR1': 'grey', 'RNR2': 'grey', 'TRNA': 'red',
               'TRNC': 'red', 'TRND':'red', 'TRNE': 'red', 'TRNF': 'red', 'TRNG': 'red',
               'TRNH': 'red', 'TRNI': 'red', 'TRNK': 'red', 'TRNL1': 'red', 'TRNL2': 'red',
               'TRNM': 'red', 'TRNN': 'red', 'TRNP': 'red', 'TRNQ': 'red', 'TRNR': 'red',
               'TRNS1': 'red', 'TRNS2':'red', 'TRNT': 'red', 'TRNV': 'red', 'TRNW': 'red',
               'TRNY': 'red'}


    
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

for i in positions:
    if len(Polym[i]) != 0:
        # find the position of the variable site
        for gene in mito_genes:
            gene_position = ''
            if i in mito_genes[gene]:
                gene_position = gene
                break
        if gene_position == '' or gene_position not in gene_colors:
            for j in polymorphism[i]:
                if j > threshold:
                    ax.scatter(i, j, edgecolor = 'black', facecolor = 'black', lw = 0, s = 5, alpha = 0.8)
        else:
            for j in polymorphism[i]:
                if j > threshold:
                    ax.scatter(i, j, edgecolor = 'black', facecolor = gene_colors[gene_position], lw = 0, s = 5, alpha = 0.8)

# restrict the x and y axis to the range of data
ax.set_xlim([0, 16570])
ax.set_ylim([0, 1])
            
# set title
if sample == 'normal':
    ax.set_title('Shared heteroplasmies in germline\n', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
elif sample == 'cancer':
    ax.set_title('Shared heteroplasmies across tumor types\n', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# set y axis label
ax.set_ylabel('Polymorphism Information Content', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels([0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000], rotation = 30, ha = 'right', size = 10, fontname = 'Helvetica', family = 'sans-serif')

plt.yticks(fontsize = 10)

# set x axis label
ax.set_xlabel('Position in mitochondrial genome', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

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
    labelbottom='off',
    direction = 'out') # labels along the bottom edge are off  
  
  
 
  
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')






