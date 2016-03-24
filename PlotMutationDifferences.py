# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:24:37 2016

@author: RJovelin
"""

# use this script to make a box plot figure with the difference between
# the number of mutations in RNA and the number of mutations in DNA in each individual
# for each tumor type

# usage: PlotMutationDifferences.py [options]
# - inputfile: file with mutational differences
# - outputfile


import os
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# get summary file with mutations
summary_file = sys.argv[1]
# get the outputfile name
outputfile = sys.argv[2]


# parse file to extract mutational differences
# open file for reading
infile = open(summary_file, 'r')
# skip header
infile.readline()
# create a dict {tumor: [mutational differences]}
Mutations = {}
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get participant
        participant = line[0]
        # get mutational differences
        mutation = int(line[1])
        tumor = line[-1]
        # populate dict
        if tumor not in Mutations:
            Mutations[tumor] = []
        Mutations[tumor].append(mutation)
infile.close()


# make a box plot figure comparing read depth betweem WGS and RNAseq for each cancer

# make a list of tumor names
tumor_names = [i for i in Mutations]
tumor_names.sort()

# combine data to a single list
all_data = []
for tumor in tumor_names:
    # append list with mutational differences
    all_data.append(Mutations[tumor])
    
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# write label for y axis
ytext = ax.set_ylabel('# RNA mutations - # DNA mutations', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
xtext = ax.set_xlabel('Tumor types', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels(tumor_names, ha = 'center', size = 10, fontname = 'Arial', family = 'sans-serif')

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# use a boxplot
bp = ax.boxplot(all_data, showmeans = True, showfliers = False, widths = 0.7, labels = tumor_names, patch_artist = True) 
    
# color WGS boxes in grey
i = 0    
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')
    if i % 2 == 0:
        # WGS data, color box 
        box.set(facecolor = [0.965, 0.745, 0.906])
    else:
        box.set(facecolor = [0.741, 0.271, 0.612])
    i += 1
        
       
# change whisker color ro black
for wk in bp['whiskers']:
    wk.set(color = 'black', linestyle = '-')
    
# change color of the caps
for cap in bp['caps']:
    cap.set(color = 'black')
        
# change the color and line width of the medians
for median in bp['medians']:
    median.set(color = 'black')
        
# change the mean marker and marker
for mean in bp['means']:
    mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 4)
    

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  


# create a list with range of x-axis values
xvals = [i + 0.5 for i in range(len(tumor_names) + 1)]
# Set a buffer around the edge of the x-axis
plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])

# do not show ticks
plt.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10)
      

# do not show ticks
plt.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are on 
    colors = 'black',
    labelsize = 10)


# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
    

# add jittered data points to box plots (when data points are not too numerous)

i = 0
for j in range(len(all_data)):
    k = np.random.uniform(-0.2, 0.2, len(all_data[j])) + (i+1)
    if (i+1) % 2 == 0:
        plt.plot(k, all_data[j], 'o', markersize = 2, markeredgecolor = [0.569, 0.251, 0.659],
                 markeredgewidth = 1.5, markerfacecolor = [0.569, 0.251, 0.659], alpha = 0.4)
    else:
        plt.plot(k, all_data[j], 'o', markersize = 2, markeredgecolor = [0.698, 0.396, 0.784],
                 markeredgewidth = 1.5, markerfacecolor = [0.698, 0.396, 0.784] , alpha = 0.4)
    i += 1


# add title
plt.title('Difference between the number of\nmutations in RNA and in DNA\n', size = 10, fontname = 'Arial')  
  
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    
    
    
    
    
    
   