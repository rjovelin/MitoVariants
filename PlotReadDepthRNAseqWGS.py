# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 15:00:13 2016

@author: RJovelin
"""

# use this script to make a box plot figure comparing read depth betweem
# WGS and RNAseq for each cancer


# usage: PlotReadDepthRNAseqWGS.py [options]
# - [mean/median]: plot mean or median read depth
# - [PerIndividual/PerPosition]: plot read depth depth per individual or per position
# - [normal/tumor]: plot read depth for normal or for tumor samples

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
# get data type
DataType = sys.argv[2]
# get sample type
SampleType = sys.argv[3]


# place this script in the folder with the heteroplasmy summary files

# create a dict {participant: tumor}
TumorID = {}
if SampleType == 'tumor':
    PrepFiles = [i for i in os.listdir('../') if ('WGS' in i or 'RNAseq' in i) and 'TPOnly' in i]
elif SampleType == 'normal':
    PrepFiles = [i for i in os.listdir('../') if ('WGS' in i or 'RNAseq' in i) and 'NTOnly' in i]    
# loop over files
for filename in PrepFiles:
    # extract cancer name from file
    cancer = filename[:filename.index('_')]
    # open file reading, loop over file, get participant IDs
    infile = open('../' + filename)
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            participant = line[0]
            TumorID[participant] = cancer
    infile.close()
print('matched participant and tumor type')

# create a dict of dict to store read depth per position for each sample
# {participant: {'RNAseq': [read depth, read depth...]}, {'WGS': [read depth, read depth...]}}
ReadDepth = {}

# make a list of tumor types
TumorTypes = ['COAD', 'LGG', 'LIRI', 'OV', 'CESC', 'SARC', 'STAD', 'UCEC', 'RECA']

# loop over cancer
for folder in TumorTypes:
    print(folder)
    # get tumor name
    tumor_name = folder
    # make a list with rnaseq or wgs subfolders
    if SampleType == 'normal':
        subdir = ['NT_RNASeq', 'NT_WGS']    
    elif SampleType == 'tumor':
        subdir = ['TP_RNASeq', 'TP_WGS']
    # loop over subdir:
    for directory in subdir:
        # create a list of subfolders with mitoseek outputs
        if 'RNASeq' in directory:
            subfolders = [i for i in os.listdir('../' + folder + '/' + directory) if 'RNASEQ' in i]
        elif 'WGS' in directory:
            subfolders = [i for i in os.listdir('../' + folder + '/' + directory) if 'GRCH37' in i]
        # removes files if files are inluded
        to_delete = []
        for i in subfolders:
            try:
                os.listdir('../' + folder + '/' + directory + '/' + i)
            except:
                to_delete.append(i)
        if len(to_delete) != 0:
            for i in to_delete:
                subfolders.remove(i)
        print('# subfolders', directory, len(subfolders))            
        # loop over subfolders
        for i in subfolders:
            # get the participant ID
            participant = i[:i.index('_')]
            # check that mitoseek ran correctly and that output has base call file
            if 'mito1_basecall.txt' in os.listdir('../' + folder + '/' + directory + '/' + i + '/'):
                # get the position-read depth for that participant {position: number of reads}
                reads = GetReadDepth('../' + folder + '/' + directory + '/' + i + '/mito1_basecall.txt')
                #  make a list with the read counts
                ReadCounts = [reads[j] for j in reads]
                # populate dicts
                if 'RNASeq' in directory:
                    assert 'RNASEQ' in i, 'subfolder should have RNAseq data'
                    # check if participant in dict
                    if participant in ReadDepth:
                        ReadDepth[participant]['RNAseq'] = list(ReadCounts)
                    else:
                        ReadDepth[participant] = {}
                        ReadDepth[participant]['RNAseq'] = list(ReadCounts)
                elif 'WGS' in directory:
                    assert 'GRCH37' in i, 'subfolder should have WGS data'
                    # check if participant in dict
                    if participant in ReadDepth:
                        ReadDepth[participant]['WGS'] = list(ReadCounts)
                    else:
                        ReadDepth[participant] = {}
                        ReadDepth[participant]['WGS'] = list(ReadCounts)
            
print('got read depth for each individual')
print('including unique samples', len(ReadDepth))
# remove individuals that do not have paired samples
to_remove = [i for i in ReadDepth if len(ReadDepth[i]) != 2]
for i in to_remove:
    del ReadDepth[i]
print('removed unique samples')
print('including paired samples', len(ReadDepth))


# create a dict {tumor: [[wgs values], [rnaseq values]]}
Coverage = {}
for participant in ReadDepth:
    # get the tumor name
    tumor_name = TumorID[participant]
    # check if tumor name is key in dict
    if tumor_name not in Coverage:
        # initialize list value
        Coverage[tumor_name] = [[], []]
    # populate list values
    # check if record read depth per individual or per position
    if DataType == 'PerIndividual':
        # check if record mean or median read depth per individual
        if StatsComp == 'median':
            # add median read depth for WGS and RNAseq
            Coverage[tumor_name][0].append(np.median(ReadDepth[participant]['WGS']))
            Coverage[tumor_name][1].append(np.median(ReadDepth[participant]['RNAseq']))
        elif StatsComp == 'mean':
            Coverage[tumor_name][0].append(np.mean(ReadDepth[participant]['WGS']))
            Coverage[tumor_name][1].append(np.mean(ReadDepth[participant]['RNAseq']))    
    elif DataType == 'PerPosition':
        Coverage[tumor_name][0].extend(list(ReadDepth[participant]['WGS']))    
        Coverage[tumor_name][1].extend(list(ReadDepth[participant]['RNAseq']))
print('recorded read depth per tumor', len(Coverage))
for tumor_name in Coverage:
    print(tumor_name, len(Coverage[tumor_name][0]), len(Coverage[tumor_name][1]))   
                

# make a box plot figure comparing read depth betweem WGS and RNAseq for each cancer

# make a list of tumor names
TumorNames = [i for i in Coverage]
TumorNames.sort()

# combine data to a single list
all_data = []
for name in TumorNames:
    # append list with wgs values
    all_data.append(Coverage[name][0])
    # append list with RNAseq values
    all_data.append(Coverage[name][1])

# create figure
fig = plt.figure(1, figsize = (11, 5))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# use a boxplot
bp = ax.boxplot(all_data, showmeans = True, showfliers = False, widths = 0.5, positions = [0, 0.5,
                                                                                           1.2, 1.7,
                                                                                           2.4, 2.9,
                                                                                           3.6, 4.1,
                                                                                           4.8, 5.3,
                                                                                           6, 6.5,
                                                                                           7.2, 7.7,
                                                                                           8.4, 8.9,
                                                                                           9.6, 10.1], patch_artist = True) 
    
# color WGS boxes in grey
i = 0    
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')
    if i % 2 == 0:
        # WGS data, color box 
        box.set(facecolor = '#1f78b4')
    else:
        box.set(facecolor = '#b2df8a')
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


# change font name

# write label for y axis
ytext = ax.set_ylabel('Read depth', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
xtext = ax.set_xlabel('Tumor types', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# create a list of tick positions
xtickpos = [0.5, 1.7, 2.9, 4.1, 5.3, 6.5, 7.7, 8.9, 10.1]

## add labels to x-ticks, rotate and align right, set size to 14
#ax.set_xticklabels(names, rotation = 30, ha = 'right', size = 10, fontname = 'Arial', family = 'sans-serif')

plt.xticks(xtickpos, TumorNames, fontsize = 10, fontname = 'Arial')


# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))


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
if DataType == 'PerIndividual':
    i = 0
    for j in range(len(all_data)):
        k = np.random.uniform(-0.2, 0.2, len(all_data[j])) + (i+1)
        if (i+1) % 2 == 0:
            plt.plot(k, all_data[j], 'o', markersize = 2, markeredgecolor = '#1f78b4',
                     markeredgewidth = 1.5, markerfacecolor = '#1f78b4', alpha = 0.4)
        else:
            plt.plot(k, all_data[j], 'o', markersize = 2, markeredgecolor = '#b2df8a',
                     markeredgewidth = 1.5, markerfacecolor = '#b2df8a', alpha = 0.4)
        i += 1


## annotate figure to add significance
## get the x and y coordinates
#yvalues = []
#for i in range(len(all_data)):
#    yvalues.extend(all_data[i])
## get max and min y values
#y_min, y_max = min(yvalues), max(yvalues)
#
## compare read depth between WGS and RNAseq for each cancer
#for i in range(0, len(all_data), 2):
#    # get the P value of Wilcoxon rank sum test
#    Pval = stats.ranksums(all_data[i], all_data[i+1])[1]
#    print(i, i+1, Pval)
#    # get stars for significance
#    if Pval > 0.05:
#        P = 'N.S.'
#    elif Pval < 0.05 and Pval > 0.01:
#        P = '*'
#    elif Pval < 0.01 and Pval > 0.001:
#        P = '**'
#    elif Pval < 0.001:
#        P = '***'
#    
#    # add bracket
#    ax.annotate("", xy=(i+1, y_max + 200), xycoords='data', xytext=(i+2, y_max + 200),
#                textcoords='data', arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
#                connectionstyle="bar,fraction=0.2", linewidth = 0.75))
#    # add stars for significance
#    if P == 'N.S.':
#        ax.text(i + 1.5, y_max + abs(y_max - y_min)*0.08, P, horizontalalignment='center',
#                verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)
#    else:
#        ax.text(i + 1.5, y_max + abs(y_max - y_min)*0.05, P, horizontalalignment='center',
#                verticalalignment='center', color = 'grey', fontname = 'Arial')



# add title
plt.title('Read depth between WGS and RNAseq\n', size = 10, fontname = 'Arial')  
  
# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')
    
    
    
    
    
    
   