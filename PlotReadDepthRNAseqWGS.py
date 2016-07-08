# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:00:13 2016

@author: RJovelin
"""

# use this script to make a box plot figure comparing read depth betweem
# WGS and RNAseq for each cancer


# usage: PlotReadDepthRNAseqWGS.py [options]
# - [mean/median]: plot mean or median read depth
# -[PerIndividual/PerPosition]: plot read depth depth per individual or per position

# import matplotlib and change api to use on server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import built in modules
import os
import sys
import numpy as np
from scipy import stats
# import custom modules
from mito_mutations import *


# get stats comparison
StatsComp = sys.argv[1]
# get date type
DataType = sys.argv[2]

# place this script in the folder with the heteroplasmy summary files


# create a dict of dict to store mean and median read depth per individual and sample type
# {participant: {'RNAseq': [mean_read_depth, median_read_depth, tumor]}, {'WGS': [mean_read_depth, median_read_depth, tumor]}}
ReadDepth = {}

# make a list of tumor types
TumorTypes = ['COAD', 'LGG', 'LIRI', 'OV', 'CESC', 'SARC', 'STAD', 'UCEC', 'RECA']

# loop over cancer
for folder in TumorTypes:
    print(folder)
    # get tumor name
    tumor_name = folder
    # create a list of subfolders
    subfolders = [i for i in os.listdir('../' + folder) if 'GRCH37' in i or 'RNASEQ' in i]
    # removes files if files are inluded
    to_delete = []
    for i in subfolders:
        try:
            os.listdir('../' + folder + '/' + i)
        except:
            to_delete.append(i)
    if len(to_delete) != 0:
        for i in to_delete:
            subfolders.remove(i)
    print('# subfolders', len(subfolders))            
    # loop over subfolders
    for i in subfolders:
        # get the participant ID
        participant = i[:i.index('_')]
        # get the position-read depth for that participant
        reads = GetReadDepth('../' + folder + '/' + i + '/mito1_basecall.txt')
        #  get the read counts
        ReadCounts = [reads[j] for j in reads]
        # compute mean and median read counts
        mean_reads = np.mean(ReadCounts)
        median_reads = np.median(ReadCounts)
        # populate dicts
        if 'GRCH37' in i:
            # check if participant in dict
            if participant in ReadDepth:
                ReadDepth[participant]['WGS'] = [mean_reads, median_reads, tumor_name]
            else:
                ReadDepth[participant] = {}
                ReadDepth[participant]['WGS'] = [mean_reads, median_reads, tumor_name]
        elif 'RNASEQ' in i:
            # check if participant in dict
            if participant in ReadDepth:
                ReadDepth[participant]['RNAseq'] = [mean_reads, median_reads, tumor_name]
            else:
                ReadDepth[participant] = {}
                ReadDepth[participant]['RNAseq'] = [mean_reads, median_reads, tumor_name]
        
print('computed mean and median read depth')
print('including unique samples', len(ReadDepth))
# remove individuals that do not have paired samples
to_remove = [i for i in ReadDepth if len(ReadDepth) != 2]
for i in to_remove:
    del ReadDepth[i]
print('removed unique samples')
print('including paired samples', len(ReadDepth))


# create a dict {tumor: [[wgs values], [rnaseq values]]}
Coverage = {}
for participant in ReadDepth:
    # get the tumor name
    tumor_name = ReadDepth[participant]['WGS'][-1]
    assert tumor_name == ReadDepth[participant]['RNAseq'][-1], 'tumor names do not match'
    # check if tumor name in dict
    if tumor_name in Coverage:
        # check if plot median or mean
        








ReadDepth = {}
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get participant
        participant = line[0]
        # check data type
        if data_type == 'mean':
            # record mean read depth
            data = float(line[1])
        elif data_type == 'median':
            data = float(line[2])
        # get tumor
        tumor = line[3]
        # get sample type
        sample_type = line[-1]
        # set up boolean
        RecordParticipant = False
        # check if using all individuals or individuals with paired data
        if individuals == 'paired':
            if participant in PairedSamples[tumor]:
                RecordParticipant = True
        elif individuals == 'all':
            RecordParticipant = True
        # check if participant is to be recorded
        if RecordParticipant:
            # populate dict
            # initialize key-value pairs 
            if tumor not in ReadDepth:
                ReadDepth[tumor] = [[], []]
            # check sample_type
            if sample_type == 'WGS':
                ReadDepth[tumor][0].append(data)
            elif sample_type == 'RNAseq':
                ReadDepth[tumor][1].append(data)
infile.close()            
            
        
# make a box plot figure comparing read depth betweem WGS and RNAseq for each cancer


# make a list of tumor names
tumor_names = [i for i in ReadDepth]
tumor_names.sort()

# combine data to a single list
all_data = []
for tumor in tumor_names:
    # append list with wgs values
    all_data.append(ReadDepth[tumor][0])
    # append list with RNAseq values
    all_data.append(ReadDepth[tumor][1])

# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# write label for y axis
ytext = ax.set_ylabel('Read depth', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
xtext = ax.set_xlabel('Tumor types', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

names = []
for tumor in tumor_names:
    names.append(tumor + '_WGS')
    names.append(tumor + '_RNA')


# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels(names, rotation = 30, ha = 'right', size = 10, fontname = 'Arial', family = 'sans-serif')


# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


# use a boxplot
bp = ax.boxplot(all_data, showmeans = True, showfliers = False, widths = 0.7, labels = names, patch_artist = True) 
    
# color WGS boxes in grey
i = 0    
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')
    if i % 2 == 0:
        # WGS data, color box 
        box.set(facecolor = [1, 0.722, 0.384])
    else:
        box.set(facecolor = [0.753, 0.949, 0.467])
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
    labelsize = 10)


# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
    

# add jittered data points to box plots (when data points are not too numerous)

i = 0
for j in range(len(all_data)):
    k = np.random.uniform(-0.2, 0.2, len(all_data[j])) + (i+1)
    if (i+1) % 2 == 0:
        plt.plot(k, all_data[j], 'o', markersize = 2, markeredgecolor = [0, 0.518, 0.176],
                 markeredgewidth = 1.5, markerfacecolor = [0, 0.518, 0.176], alpha = 0.4)
    else:
        plt.plot(k, all_data[j], 'o', markersize = 2, markeredgecolor = [0.957, 0.631, 0],
                 markeredgewidth = 1.5, markerfacecolor = [0.957, 0.631, 0] , alpha = 0.4)
    i += 1


# annotate figure to add significance
# get the x and y coordinates
yvalues = []
for i in range(len(all_data)):
    yvalues.extend(all_data[i])
# get max and min y values
y_min, y_max = min(yvalues), max(yvalues)

# compare read depth between WGS and RNAseq for each cancer
for i in range(0, len(all_data), 2):
    # get the P value of Wilcoxon rank sum test
    Pval = stats.ranksums(all_data[i], all_data[i+1])[1]
    # get stars for significance
    if Pval > 0.05:
        P = 'N.S.'
    elif Pval < 0.05 and Pval > 0.01:
        P = '*'
    elif Pval < 0.01 and Pval > 0.001:
        P = '**'
    elif Pval < 0.001:
        P = '***'
    
    # add bracket
    ax.annotate("", xy=(i+1, y_max + 200), xycoords='data', xytext=(i+2, y_max + 200),
                textcoords='data', arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                connectionstyle="bar,fraction=0.2", linewidth = 0.75))
    # add stars for significance
    if P == 'N.S.':
        ax.text(i + 1.5, y_max + abs(y_max - y_min)*0.08, P, horizontalalignment='center',
                verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)
    else:
        ax.text(i + 1.5, y_max + abs(y_max - y_min)*0.05, P, horizontalalignment='center',
                verticalalignment='center', color = 'grey', fontname = 'Arial')



# add title
plt.title('Read depth between WGS and RNAseq\n', size = 10, fontname = 'Arial')  
  
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    
    
    
    
    
    
   