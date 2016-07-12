# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 11:54:09 2016

@author: RJovelin
"""

# use this script to plot correlation between genotype at MRPP3 and RNA heteroplasmy frequency

# usage PlotCorrelMRPP3GenotypeAltFreq [options]
# -[trna/P9/allgenes]: consider only RNA modifications in tRNA, modifs in tRNA P9 or all genes

# place this script in folder HeteroPlasmyFilesRNAOnly with summary files

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
from scipy import stats
# import custom modules
from mito_mutations import *
print('done importing modules')

# fix threshold to detect heteroplasmy to 1%
threshold = 1

# get site type from command
SiteType = sys.argv[1]


# create a dict {individual: genotype}
Genotypes = {}

# make a list of tumor types with VCF files
VCFTumor = ['COAD', 'OV', 'RECA-EU', 'UCEC', 'CESC', 'LGG', 'LIRI', 'SARC', 'STAD']

# loop over tumor folders, unzip VCF files
for folder in VCFTumor:
    print(folder)
    # make a list of gzip files    
    if folder == 'STAD':
        files = [i for i in os.listdir('../../../VCFs/Germlines_vcf/perCancerType/' + folder + '/filteredVCFS') if i[-3:] == '.gz']
    else:
        files = [i for i in os.listdir('../../../VCFs/Germlines_vcf/perCancerType/' + folder) if i[-3:] == '.gz']
    if len(files) != 0:
        # unzip files
        if folder == 'STAD':
            for filename in files:
                os.system('gunzip ' + '../../../VCFs/Germlines_vcf/perCancerType/' + folder + '/filteredVCFS/' + filename)
        else:
            for filename in files:
                os.system('gunzip ' + '../../../VCFs/Germlines_vcf/perCancerType/' + folder + '/' + filename)
print('done unzipping VCF files')                

# make a list of WGS normal files with corresponding individual ID and 
MatchFiles = [i for i in os.listdir('../') if 'WGS_NTOnly' in i]
assert len(MatchFiles) == 9, 'there should be matching files for 9 tumors'

# create a dict {BAM ID : individual ID} 
MatchingIDs = {}
for filename in MatchFiles:
    print(filename)
    infile = open('../' + filename)
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            # get individual ID
            Individual = line[0]
            # get BAM ID
            BamID = line[2]
            BamID = BamID[BamID.index('wgsmito/') + len('wgsmito/'): BamID.index('.mito.')]
            # populate dict
            MatchingIDs[BamID] = Individual
    infile.close()
print('matched Bams with Individual IDs', len(MatchingIDs))

# loop over folders, extract MRPP3 genotypes of all individuals
for folder in VCFTumor:
    print(folder)
    if folder == 'STAD':
         mrpp3 = GetMRPP3Genotypes('../../../VCFs/Germlines_vcf/perCancerType/' + folder + '/filteredVCFS/')
    else:
        mrpp3 = GetMRPP3Genotypes('../../../VCFs/Germlines_vcf/perCancerType/' + folder + '/')
    # populate Genotypes dict with individual ID
    for ID in mrpp3:
        if ID in MatchingIDs:
            Genotypes[MatchingIDs[ID]] = mrpp3[ID]
print('extracted MRPP3 genotype', len(Genotypes))                

# get the heteroplasmy frequency for trna, p9 or all genes 

# get the positions of each mitochiondrial gene
mito_genes = MitoAnnotation('rCRS_genes_MT.text.txt')
# make a set of tRNA positions
trna_indices = []
for gene in mito_genes:
    if gene.startswith('TRN'):
        for i in mito_genes[gene]:
            trna_indices.append(i)
print('got tRNA indices')

# get the gene coordinates [start, end, orientation]
mito_coords = MitoCoordinates('rCRS_genes_MT.text.txt')
print('got gene coordinates')


# make a list of summary files 
SummaryFiles = [i for i in os.listdir() if 'tumor_RNAOnly' in i and '.txt' in i]

# create a dict {mutation: [list of allele frequencies]}
mutations = {}

# loop over filename in files
for filename in files:
    # open file for reading
    infile = open(filename)
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # grab participant ID
            participant = line[1]
            # get position 0-based
            position = int(line[0]) - 1
            # get ref allele
            ref = line[6]
            # get major and minor allele
            major, minor = line[7], line[8]
            # get major and minor counts
            major_counts, minor_counts = int(line[9]), int(line[10])
            # get gene
            gene = line[2]
            # verify that major and minor read counts > threshold
            alleles = [major, minor]
            total_reads = 0
            allele_freq = []
            for i in range(12, 20):
                total_reads += int(line[i])
            for i in range(len(alleles)):
                if alleles[i] == 'A':
                    allele_counts = int(line[12]) + int(line[16])
                elif alleles[i] == 'T':
                    allele_counts = int(line[13]) + int(line[17])
                elif alleles[i] == 'C':
                    allele_counts = int(line[14]) + int(line[18])
                elif alleles[i] == 'G':
                    allele_counts = int(line[15]) + int(line[19])
                assert (allele_counts / total_reads) * 100 > threshold, 'allele frequency is lower than heteroplasmy threshold'   
                # check if allele counts is major or minor allele count                
                if i == 0:
                    assert allele_counts == major_counts, 'allele count is different than major count'
                elif i == 1:
                    assert allele_counts == minor_counts, 'allele count is different than minor count'
                allele_freq.append(allele_counts / total_reads)
            # identify mutant allele by comparing minor and major to ref
            # get the frequency (read count / total reads at position) for the mutant allele
            if major == ref and minor != ref:
                # mutant is minor
                freq = allele_freq[1]
            elif major != ref and minor == ref:
                # mutant is major            
                freq = allele_freq[0]
            # check which sites need to be recorded
            if SiteType == 'allgenes':
                # record all sites
                if participant in mutations:
                    mutations[participant].append(freq)
                else:
                    mutations[participant] = [freq]
            elif SiteType == 'trna':
                # record tRNA sites
                if line[2].startswith('TRN'):
                    assert position in trna_indices, 'position should be trna'
                    # position is in tRNA, populate dict
                    if participant in mutations:
                        mutations[participant].append(freq)
                    else:
                        mutations[participant]  [freq]
            elif SiteType == 'P9':
                assert position in trna_indices, 'position should be a trna index'
                # record tRNA P9 sites
                # find the tRNA corresponding to current site
                for i in mito_genes:
                    if position in mito_genes[i]:
                        # stop looking, found tRNA gene
                        break
                assert i == gene, 'gene name does not match woth expected position'
                # convert genomic position to tRNA position
                assert position in range(mito_coords[gene][0], mito_coords[gene][1]), 'position should be in gene'
                trnaposition = GenomicPositionToGenePosition(position, mito_coords[gene][0], mito_coords[gene][1], mito_coords[gene][2])
                # record only P9 posiitons (8 in 0-based index)
                if trnaposition == 8:
                    if participant in mutations:
                        mutations[participant].append(freq)
                    else:
                        mutations[participant] = [freq]
                 
    infile.close()
print('got allele frequencies')            

# create lists with frequencies for each geneotype
AA = [freq for freq in mutations[participant] if Genotypes[participant] == 'AA']
AG = [freq for freq in mutations[participant] if Genotypes[participant] == 'AG']
GG = [freq for freq in mutations[participant] if Genotypes[participant] == 'GG']        
        
# count the number of participants
total_participants = 0
# loop over individual IDs in frequency dict
for participant in mutations:
    # check if participant has genotype
    if participant in Genotypes:
        # update counter
        total_participants += 1
print('total participants', total_participants)
print('AA', len(AA))
print('AG', len(AG))
print('GG', len(GG))


# create a figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)

# make a list with all data
Data = [AA, AG, GG] 

# plot jittered data points 
i = 0
for j in range(len(Data)):
    k = np.random.uniform(-0.2, 0.2, len(Data[j])) + (i+1)
    plt.plot(k, Data[j], 'o', markersize = 2, markeredgecolor = 'black', markeredgewidth = 1.5, markerfacecolor = 'black', alpha = 0.5)
    i += 1

# limit x and y axis
ax.set_ylim([0, 1])
ax.set_xlim([0,4])
  
# set font for all text in figure
FigFont = {'fontname':'Helvetica'}

# set y axis label
ax.set_ylabel('Alternative Allele Frequency', size = 10, ha = 'center', **FigFont)

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels(['A/A', 'A/G', 'G/G'], ha = 'center', size = 10, **FigFont)

# set x axis label
ax.set_xlabel('rs11156878', size = 10, ha = 'center', **FigFont)

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
    label.set_fontname('Helvetica')
    

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

