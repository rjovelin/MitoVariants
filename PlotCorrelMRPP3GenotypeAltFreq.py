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
    # make a list of gzip files    
    if folder == 'STAD':
         mrpp3 = GetMRPP3Genotypes('../../../VCFs/Germlines_vcf/perCancerType/' + folder + '/filteredVCFS')
    else:
        mrpp3 = GetMRPP3Genotypes('../../../VCFs/Germlines_vcf/perCancerType/' + folder)
    # populate Genotypes dict with individual ID
    ID = [i for i in mrpp3]
    assert len(ID) == 1, 'there should be a single indivual'
    ID = ID[0]
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










# sort frequency values
for effect in mutations:
    mutations[effect].sort()

# make list with functional category and data
Data = [[key, val] for key, val in mutations.items()]
# replace some functional categories with a shorter name
for i in range(len(Data)):
    if Data[i][0] == 'stopgain':
        Data[i][0] = 'PSC'
    elif Data[i][0] == 'stoploss':
        Data[i][0] = 'SCL'
    elif Data[i][0] == 'non-synonymous':
        Data[i][0] = 'NonSyn'
    elif Data[i][0] == 'synonymous':
        Data[i][0] = 'Syn'
    # keep the same name for the other categories (tRNA, DLoop, NonCoding, Ribosomal)

# sort according to names
Data.sort()

# create parallel lists with functional names and data
categories = [Data[i][0] for i in range(len(Data))]
data = [Data[i][1] for i in range(len(Data))]    
print(categories)


# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# create a list to store the variables for building the legend
Graphs = []

# make a list of colors
colorscheme = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5']

# loop over categories
for i in range(len(categories)):
    print(categories[i])
    graph = ax.step(data[i], np.linspace(0, 1, len(data[i]), endpoint=True), linewidth = 1.2, color = colorscheme[i], alpha = 0.7)
    Graphs.append(graph)
print('plotted CDF')

# add label for the Y axis
ax.set_ylabel('Cumulative fraction of mutations', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Mutant allele Frequency', size = 10, ha = 'center', fontname = 'Arial')

# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# do not show ticks on 1st graph
ax.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

# do not show ticks
ax.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# add lines
lns = Graphs[0]
for i in range(1, len(Graphs)):
    lns = lns + Graphs[i]
# get labels
labs = [categories[i] for i in range(len(categories))]
print(labs)
# plot legend
ax.legend(lns, labs, loc=4, fontsize = 6, frameon = False)

# build outputfile namewith parameters
# extract the cancer name
if which_files == 'singlefile' and sample == 'tumor':
    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_' + sample)]
elif which_files == 'singlefile' and sample == 'specific':
    cancer = HeteroSummaryFile[HeteroSummaryFile.index('_') + 1: HeteroSummaryFile.index('_TumorSpecific')]

if which_files == 'singlefile':
    outputfile = 'CDFMutantAlleleFreq' + cancer + sample.capitalize() + '.pdf'
elif which_files == 'allfiles':
    outputfile = 'CDFMutantAllelefreq' + sample.capitalize() + '.pdf'

fig.savefig(outputfile, bbox_inches = 'tight')







#####################################




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




