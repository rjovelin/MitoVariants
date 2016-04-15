# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 15:12:35 2016

@author: RJovelin
"""

# use this script to generate a summary file with sample frequencies of multiallelic sites
# for each tumor type

# usage: python3 MakeSummaryMultiAllelicFile.py [options]
# - HeteroplasmyFile: summary file with heteroplasmies for a given sample
# - rCRS: mitochondrial genome annotation
# - threshold: percent of reads supporting alternative alleles
# - outputfile

import sys
from mito_mutations import *

# get heteroplasmy file
HeteroplasmyFile = sys.argv[1]
# get mitochindrial annotation
rCRS = sys.argv[2]
# get threshold
threshold = float(sys.argv[3])
# get outputfile name
outputfile = sys.argv[4]

# extract the tumor name from the Heteroplasmy file name
cancer = HeteroplasmyFile[HeteroplasmyFile.index('_') + 1 : HeteroplasmyFile.index('_', HeteroplasmyFile.index('_') + 1)] 

# check if file already exists
try:
    newfile = open(outputfile, 'r')
except:
    # write header if file does not exist
    newfile = open(outputfile, 'w')
    newfile.write('\t'.join(['tumor', 'position', 'N_alleles', 'alleles', 'gene', 'position_in_gene', 'N', 'sample_size', 'sample_frequency']) + '\n')
    print('content written to new file')
else:
    # if it exists, append content
    newfile = open(outputfile, 'a')
    print('content appended to existing file')

# identify multiallelic SNPs
# get the number of alleles at each position for each individual 
# {individual : {position : [allele1, allele2, allele3, allele4]}}
AlleleCounts = IdentifyMultiAllelicPositions(HeteroplasmyFile, threshold)

# get the positions of each mitochiondrial gene {gene: set(positions)}
mito_genes = MitoAnnotation(rCRS)
# get the gene coordinates 0-based
for i in mito_genes:
    mito_genes[i] = list(mito_genes[i])
    # sort positions
    mito_genes[i].sort()
    # get the start and end position 0-based {gene:[start:end]}
    # add 1 to end to convert to 0-based because end is inclusive in gene positions
    mito_genes[i] = [mito_genes[i][0], mito_genes[i][-1] + 1]

# create a dict {position: [gene, orientation, sample_size]}
annotation = {}
infile = open(HeteroplasmyFile, 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get position 0-based
        position = int(line[0]) - 1
        gene = line[2]
        sample_size = int(line[11])
        orientation = line[3]
        if '_' in orientation:
            orientation = orientation[orientation.index('_')+2: -1]
        if position not in annotation:
            annotation[position] = [gene, orientation, int(line[11])]
        else:
            assert annotation[position] == [gene, orientation, int(line[11])], 'gene, orientation or sample size don\'t match'
infile.close()

# convert dict into {position: alleles: N_individuals}
Positions = {}
# loop over individuals
for individual in AlleleCounts:
    # loop over positions
    for position in AlleleCounts[individual]:
        alleles = AlleleCounts[individual][position]
        # only consider multiallelic snps (> 2 alleles)
        if len(alleles) > 2:
            # sort alleles
            alleles.sort()
            # get alleles
            alleles = ''.join(alleles)
            # populate dict
            if position not in Positions:
                Positions[position] = {}
            if alleles in Positions[position]:
                Positions[position][alleles] += 1
            else:
                Positions[position][alleles] = 1


# create a sorted list of positions
SortedPositions = [i for i in Positions]
SortedPositions.sort()

# loop over sorted positions
for position in SortedPositions:
    # get gene coordinates and sample size
    gene, orientation, sample_size = annotation[position][0], annotation[position][1], annotation[position][2]
    if gene != 'NA':
        # check if position is found in 2 overlapping genes
        if '|' in gene:
            # grap first gene
            gene = gene[:gene.index('|')]
        start, end = mito_genes[gene][0], mito_genes[gene][1] 
        # get the position of the SNP on the gene (in 5' to 3' orientation)
        SNP_position = GenomicPositionToGenePosition(position, start, end, orientation)
        # loop over alleles
        for alleles in Positions[position]:
            # write content to file, convert position 1-based
            newfile.write('\t'.join([cancer, str(position + 1), str(len(alleles)),
                                     alleles, gene, str(SNP_position + 1), str(Positions[position][alleles]),
                                     str(sample_size), str(round((Positions[position][alleles] / sample_size)*100, 2))]) + '\n')
    elif gene == 'NA':
        # SNP in non-genic position, assign positions to NA
        start, end, SNP_position = 'NA', 'NA', 'NA' 
        # loop over alleles
        for alleles in Positions[position]:
            # write content to file, convert position 1-based
            newfile.write('\t'.join([cancer, str(position + 1), str(len(alleles)),
                                     alleles, gene, SNP_position, str(Positions[position][alleles]),
                                     str(sample_size), str(round((Positions[position][alleles] / sample_size)*100, 2))]) + '\n')

newfile.close()