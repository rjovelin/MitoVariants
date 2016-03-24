# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 15:59:35 2016

@author: RJovelin
"""

# use this script to make a box plot figure comparing read depth betweem
# heteroplasmic and invariant positions either for WGS or RNAseq for a given cancer


# usage: PlotreadDepthVariablePositions.py [options]
# - inputfile: json file with dict {participant: [[read_depth_snps], [read_depth_not_variable]]}
# - [WGS/RNAseq]: datatype
# - outputfile: figure file


import os
import sys
import numpy as np
from scipy import stats
import json
import matplotlib.pyplot as plt

# get json file with read depth
datafile = sys.argv[1]
# get data type
data_type = sys.argv[2]
# get the outputfile name
outputfile = sys.argv[3]

# get the tumor type from the datafile
tumor_type = datafile[:datafile.index('_')]

# open file for reading
infile = open(datafile, 'r')
# load json dict {participant: [[read_depth_snps], [read_depth_not_variable]]}
ReadDepth = json.load(infile)
# close file
infile.close()

# create  a data list [[read_depth_snps], [read_depth_not_variable]]
AllData = [[], []]
# loop over dict
for participant in ReadDepth:
    # dump read depth in each inner list
    AllData[0].extend(ReadDepth[participant][0])
    AllData[1].extend(ReadDepth[participant][1])


# make a box plot figure comparing read depth betweem variable and non-variable positions

# make a list of position types
names = ['Variable', 'Conserved']

# create figure
#fig = plt.figure(1, figsize = (4.3,2.56))
fig = plt.figure(1, figsize = (2.5, 3))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# write label for y axis
ytext = ax.set_ylabel('Read depth', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
xtext = ax.set_xlabel('MT sites', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels(names, ha = 'center', size = 10, fontname = 'Arial', family = 'sans-serif')

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# use a boxplot
bp = ax.boxplot(AllData, showmeans = True, showfliers = False, widths = 0.7, labels = names, patch_artist = True) 
    
# color WGS boxes in grey
i = 0    
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')    
    if i % 2 == 0:
        # SNPs data, color box in grey
        box.set(facecolor = [0.518, 0.439, 0.804])
    else:
        box.set(facecolor = [0.808, 0.765, 0.9])
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
    mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 5)
    

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  


# create a list with range of x-axis values
xvals = [i + 0.5 for i in range(len(names) + 1)]
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
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on'


# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
    
## note some ways to modify label font size and color of the ticks
#for label in ax.get_yticklabels():
#    label.set_color('black')
#    label.set_fontsize(10)

# add title
plt.title('Read depth in {0} for {1} data\n'.format(tumor_type, data_type), size = 12, fontname = 'Arial')  



# annotate figure to add significance

# get the P value of Wilcoxon rank sum test
Pval = stats.ranksums(AllData[0], AllData[1])[1]
# get stars for significance
if Pval > 0.05:
    P = 'N.S.'
elif Pval < 0.05 and Pval > 0.01:
    P = '*'
elif Pval < 0.01 and Pval > 0.001:
    P = '**'
elif Pval < 0.001:
    P = '***'
    
# get the x and y coordinates
yvalues = []
yvalues.extend(AllData[0])    
yvalues.extend(AllData[1])
# get max and min y values
y_min, y_max = min(yvalues), max(yvalues)

# add bracket
ax.annotate("", xy=(1, y_max), xycoords='data', xytext=(2, y_max),
            textcoords='data', arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
            connectionstyle="bar,fraction=0.1"))
# add stars for significance
ax.text(1.5, y_max + abs(y_max - y_min)*0.05, P, horizontalalignment='center',
       verticalalignment='center', color = 'grey', fontname = 'Arial')


# add jittered data points to box plots (when data points are not too numerous)
# set x up coordinates
#a = np.random.uniform(-0.2, 0.2, len(AllData[0])) + 1
#b = np.random.uniform(-0.2, 0.2, len(AllData[1])) + 2
#plt.plot(a, AllData[0], 'o', markersize = 4, markeredgecolor='blue', markeredgewidth = 1.5, color = 'blue', alpha = 0.3)
#plt.plot(b, AllData[1], 'o', markersize = 4, markeredgecolor='blue', markeredgewidth = 1.5, color = 'blue', alpha = 0.3)


# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    
  