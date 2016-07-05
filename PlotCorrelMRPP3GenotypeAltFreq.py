# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 11:54:09 2016

@author: RJovelin
"""

# use this script to plot correlation between genotype at MRPP3 and RNA modif frequency

# usage PlotCorrelMRPP3GenotypeAltFreq [options]
# -[altfreq, heterofreq]: plot alternative frequency in the population or heteroplasmy level within indivual


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
# import custom modules
from mito_mutations import *





# create a correlation plot between alternative allele frequency as a function of genotypes at given position

# usage plot_correlation_altfreq_genotype.py [parameters]
# inputfile : file with genotypes and allele frequencies
# outputfile : figure file
# SNP ID: ID of the SNP with genotype

import os
import sys
import matplotlib.pyplot as plt
import numpy as np




# get inputfile from command
genotype_file = sys.argv[1]
# get outputfile
figure_file = sys.argv[2]
# get the ID of the SNP
SNP_ID = sys.argv[3]

# create a dict with {genotypes: list of allelic frequencies}
genotypes = {}
# open file for reading
infile = open(genotype_file, 'r')
#loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get genotype
        genotype = line[1]
        # make a list of allelic frequencies
        freq = line[2].split(';')
        altfreq = list(map(lambda x: float(x), freq))
        # check if genotype in dict
        if genotype in genotypes:
            # add all frequency values
            genotypes[genotype].extend(altfreq)
        else:
            # add genotype : list of frequencies pairs to dict
            genotypes[genotype] = altfreq
# close file 
infile.close()


# make a list of genotypes
SNP = [i for i in genotypes]
# sort list
SNP.sort()


# create a figure
fig = plt.figure(1, figsize = (5, 3))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)

#for i in range(len(SNP)):
#    for j in genotypes[SNP[i]]:
#        ax.scatter(i+1, j, edgecolor = 'black', facecolor = 'black', lw = 1, s = 10, alpha = 0.5)



xA = np.random.normal(1, 0.1, len(genotypes[SNP[0]])) 
         

xB = np.random.normal(2, 0.1, len(genotypes[SNP[1]]))

xC = np.random.normal(3, 0.1, len(genotypes[SNP[2]]))


plt.scatter(xA, genotypes[SNP[0]], edgecolor = 'black', facecolor = 'black', lw = 1, s = 10, alpha = 0.5)
plt.scatter(xB, genotypes[SNP[1]], edgecolor = 'grey', facecolor = 'grey', lw = 1, s = 10, alpha = 0.5) 
plt.scatter(xC, genotypes[SNP[2]], edgecolor = 'lightgrey', facecolor = 'lightgrey', lw = 1, s = 10, alpha = 0.5)


ax.set_ylim([0, 1])
ax.set_xlim([0,4])
            
# set y axis label
ax.set_ylabel('Alternative Allele Frequency', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels([SNP[0], SNP[1], SNP[2]], ha = 'center', size = 10, fontname = 'Helvetica', family = 'sans-serif')

# set x axis label
ax.set_xlabel(SNP_ID, size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


# set positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+1 for i in range(3)]

plt.xticks(tick_pos, SNP)














# save figure
fig.savefig(figure_file, bbox_inches = 'tight')




