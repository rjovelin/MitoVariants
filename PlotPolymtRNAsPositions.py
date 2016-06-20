# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 12:50:38 2016

@author: RJovelin
"""

# create a figure with polymorphism content at each positions of tRNAs

# place this script in folder with heteroplasmy summary files


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


# usage python3 PlotPolymtRNAsPositions.py [parameters]
# - [singlefile/allfiles]: whether a single summary file or multiple summary files
# - [tumor/specific]: whether RNA variants are in tumor or are tumor specific (filtered based on normale)
# - threshold: the value above which PIC is considered (eg. threshold = 0 means 0 values won't be shown)
# - [True/False] : whether sites with singleton variants are removed or not 
# - HeteroplasmySummaryFile: summary file of tumor if singlefile is used


#use single of all summary files
which_files = sys.argv[1]
sample = sys.argv[2]
threshold = float(sys.argv[3])
RemoveSingleton = sys.argv[4]
if RemoveSingleton == 'True':
    RemoveSingleton = True
elif RemoveSingleton == 'False':
    RemoveSingleton = False


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

# create a dict to store the polymorphism values at each tRNA positions
tRNAPolym = {}

# loop over files
for filename in files:
    # compute polymorphism for each position
    polymorphism = GenomePositionsPic(filename)
    # remove positions that are not tRNA
    to_delete = [i for i in polymorphism if i not in trna_indices]
    for i in to_delete:
        del polymorphism[i]
    # convert genomic positions to tRNA positions
    for i in polymorphism:
        # find the tRNA for that position
        for gene in mito_genes:
            if i in mito_genes[gene]:
                # stop looking, got the gene
                break
        # check that correct gene is found
        assert i in range(mito_coords[gene][0], mito_coords[gene][1]), 'position should be in {0}'.format(gene)
        # convert coordinate to trna coordinate
        position = GenomicPositionToGenePosition(i, mito_coords[gene][0], mito_coords[gene][1], mito_coords[gene][2])
        # populate dict
        if position in tRNAPolym:
            tRNAPolym[position].extend(polymorphism[i])
        else:
            tRNAPolym[position] = list(polymorphism[i])
        
        
# create a list of positions
positions = [i for i in tRNAPolym]


    
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

for i in positions:
    # set up boolean
    RecordSite = False
    if RemoveSingleton == True:
        if len(polymorphism[i]) > 1:
            RecordSite = True
    elif RemoveSingleton == False:
        if len(polymorphism[i]) != 0:
            RecordSite = True
    if RecordSite == True:
        # find the position of the variable site
        for gene in mito_genes:
            if i in mito_genes[gene]:
                break
        if gene == '' or gene not in gene_colors:
            for j in polymorphism[i]:
                if j > threshold:
                    ax.scatter(i, j, edgecolor = 'black', facecolor = 'black', lw = 0, s = 5, alpha = 0.8)
        else:
            for j in polymorphism[i]:
                if j > threshold:
                    ax.scatter(i, j, edgecolor = 'black', facecolor = gene_colors[gene], lw = 0, s = 5, alpha = 0.8)

# restrict the x and y axis to the range of data
ax.set_xlim([0, 16570])
ax.set_ylim([0, 1])
            
# set title
if sample == 'normal':
    ax.set_title('{0} - Heteroplasmies in germline\n'.format(cancer), size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
elif sample == 'tumor':
    ax.set_title('{0} - Heteroplasmies in tumor\n'.format(cancer), size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
elif sample == 'specific':
    ax.set_title('{0} - Tumor-specific heteroplasmies\n'.format(cancer), size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

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
    labelbottom='off') # labels along the bottom edge are off  
  
# build outputfile with parameters  
if RemoveSingleton == True:
    outputfile = 'Polymorphism' + cancer + sample + 'NoSingletons' + '.pdf' 
elif RemoveSingleton == False:
    outputfile = 'Polymorphism' + cancer + sample + '.pdf' 
  
 
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')

