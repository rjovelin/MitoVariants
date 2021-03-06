# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:20:37 2016

@author: RJovelin
"""

# use this script to plot polymorphism content for all tumors types pooled together

# import matplotlib and change api to use on server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import built in modules
import sys
# import custom modules
from mito_mutations import *


# place script in same folder with heteroplasmy summary files

# usage python3 PlotPolymAllTumors.py [parameters]
# - threshold: the value above which PIC is considered (eg. threshold = 0 means 0 values won't be shown)
# - [tumor/specific]: RNA variants are in tumor or are tumor specific (filtered with normal) 
# - [full/limited/black] : color panel: each gene has a different color (full) 
#                        or proteins are in grey and tRNAs in red (limited)
#                        or only tRNAs are labeled in red (black)


# get thresold from command
threshold = float(sys.argv[1])
# get normal or cancer from command
sample = sys.argv[2]
# get the color panel to use
colors = sys.argv[3]


# make a list of hetero summary files
if sample == 'tumor':
    files = [i for i in os.listdir() if sample in i and 'RNAOnly' in i and '.txt' in i]
elif sample == 'specific':
    files = [i for i in os.listdir() if 'TumorSpecific' in i and '.txt' in i] 

# create a set of SNPs among all tumor types
AllSnps = PoolVariablePositions(files)
print('all snps', len(AllSnps))

# get the positions of each mitochiondrial genes and regions
mito_genes = MitoAnnotation('rCRS_genes_MT.text.txt')

# create a dict to with common snps : list of PIC values pairs
Polym = {}
for i in AllSnps:
    Polym[i] = []

# loop over summary file
for filename in files:
    # compute PIC for each individual at each position
    polymorphism = GenomePositionsPic(filename)
    print(filename, len(polymorphism))
    # add PIC values to dict
    for i in polymorphism:
        if i in Polym:
            Polym[i].extend(list(polymorphism[i]))
        
        
# create a list of positions
positions = [i for i in Polym]
positions.sort()

# check which color panel to use
# create a dict with gene : color
if colors == 'full':
    # tRNAs are in red, each portein coding gene has a different color    
    gene_colors = {'ATP6': '#8dd3c7', 'ATP8': '#bc80bd',
                   'COX1': '#ccebc5', 'COX2': '#fdb462', 
                   'COX3': '#b3de69', 'CYTB': '#fccde5',
                   'ND1': '#ffffb3', 'ND2': '#ffed6f',
                   'ND3': '#fb8072', 'ND4': '#80b1d3',
                   'ND4L': '#80b1d3', 'ND5': '#b3de69',
                   'ND6': '#fdb462', 'RNR1': '#d9d9d9',
                   'RNR2': '#bebada', 'TRNA': 'red',
                   'TRNC': 'red', 'TRND':'red', 'TRNE': 'red', 'TRNF': 'red', 'TRNG': 'red',
                   'TRNH': 'red', 'TRNI': 'red', 'TRNK': 'red', 'TRNL1': 'red', 'TRNL2': 'red',
                   'TRNM': 'red', 'TRNN': 'red', 'TRNP': 'red', 'TRNQ': 'red', 'TRNR': 'red',
                   'TRNS1': 'red', 'TRNS2':'red', 'TRNT': 'red', 'TRNV': 'red', 'TRNW': 'red',
                   'TRNY': 'red'}

elif colors == 'black':
    # tRNAs are in red all others are in black
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
    # tRNAs are in red, genes are in grey and noncoding are in black
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
            if i in mito_genes[gene]:
                break
        if gene == '' or gene not in gene_colors:
            for j in Polym[i]:
                if j > threshold:
                    ax.scatter(i, j, edgecolor = 'black', facecolor = 'black', lw = 0, s = 5, alpha = 0.8)
        else:
            for j in Polym[i]:
                if j > threshold:
                    ax.scatter(i, j, edgecolor = 'black', facecolor = gene_colors[gene], lw = 0, s = 5, alpha = 0.8)

# restrict the x and y axis to the range of data
ax.set_xlim([0, 16570])
ax.set_ylim([0, 1])
            
# set title
if sample == 'tumor':
    ax.set_title('Pooled heteroplasmies across tumor types\n', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
elif sample == 'specific':
    ax.set_title('Pooled tumor-specific heteroplasmies\n', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

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
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off',
    direction = 'out') # labels along the bottom edge are off  
  
  
# build outputfile with parameters
outputfile = 'PooledPolymorphism' + sample.capitalize() + '.pdf' 
  
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
