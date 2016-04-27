# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 11:16:15 2016

@author: RJovelin
"""

# import modules
import os
import numpy as np



# use this function to make a list of directories with mitoseek outputs
def FindDuplicateIndividuals(folder):
    '''
    (str) -> list
    Look for all the subdirectories in folder and return a list of
    subdirectory names to ignore because of duplicate mapping to GRCH37 and
    GRCH37lite references
    '''
    
    # create a list of subfolder, exluding any files
    subfolders = os.listdir(folder)
    to_delete = []
    for i in subfolders:
        try:
            os.listdir(i)
        except:
            to_delete.append(i)
    for i in to_delete:
        subfolders.remove(i)
    
    # create a list with individual names
    names = []
    for i in subfolders:
        # loop over folder names in folder
        if 'GRCH37' in i:
            # parse name and extract participant name
            participant = i[: i.index('_GRCH37')]
            names.append(participant)
    # create a list with folder names to ignore
    duplicates = []
    # loop over file names
    for i in names:
        # check if participant is mapped to both references
        if names.count(i) > 1:
            # check that participant ID is counted twice
            assert names.count(i) == 2, 'Participant ID is counted more than twice'
            # give preference to GRCH37 reference
            # check that duplicate ref is is folder
            assert i + '_GRCH37lite' in subfolders, 'Participant not mapped to GRCH37'
            # add GRCH37 ref to duplicate folder
            duplicates.append(i + '_GRCH37lite')
    return duplicates
    



# use this function to make a list of folders corresponding to GRCH37 ref
def GetValidReference(folder):
    '''
    (str) -> list
    Look for all the subdirectories in folder and return a list of
    subdirectory names with GRCH37 as reference
    '''
    
    # create a list of subfolder, exluding any files
    subfolders = os.listdir(folder)
    to_delete = []
    for i in subfolders:
        try:
            os.listdir(i)
        except:
            to_delete.append(i)
    for i in to_delete:
        subfolders.remove(i)
    
    # create a list of valid subfolders
    valid = []
    for i in subfolders:
        # loop over folder names in folder
        if i[-7:] == '_GRCH37':
            # keep subfolder with GRCH37 as ref
            valid.append(i)
    return valid    
    
    
# use this function to determine the sample size of each position
def SampleSize(folder, tissue_type, dataset, minimum_coverage):
    '''
    (str, str, str, int) -> dict
    Return a dict with the individuals with > minumum coverage at each position
    of the mitochondrial genome for the tumor or normal tissue type and for
    WGS or for RNA
    '''
    
    # use only individuals mapped with GRCH37 ref    
    if dataset == 'WGS':
        directories = GetValidReference(folder)
    elif dataset == 'RNA':
        # all individuals have been mapped using GRCH37
        # create a list of subdirectories
        directories = os.listdir()
        to_delete = []
        for i in directories:
            try:
                os.listdir(i)
            except:
                to_delete.append(i)
        for i in to_delete:
            directories.remove(i)    
    
    # create a dict with position 0-based {position: set(indidivuals)}    
    sample = {}
    
    # loop over directories
    for i in directories:
        # get participant ID
        participant = i[:i.index('_')]
        # check if mitoseek was run on paired samples or on individual samples
        if 'mito2_basecall.txt' in os.listdir(i):
            # check if tissue type is normal or tumor
            if tissue_type == 'normal':
                file = i + '/' + 'mito2_basecall.txt'
            elif tissue_type == 'tumor':
                file = i + '/' + 'mito1_basecall.txt'
        else:
            # mitoseek was ran on a single sample
            # normal or tumor is already specified by the folder
            file = i + '/' + 'mito1_basecall.txt'
        
        # open file for reading
        infile = open(file, 'r')
        # skip header
        infile.readline()
        # loop over file, grab position, update with count
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split('\t')
                position = int(line[1]) - 1
                # get read depth
                coverage = sum(list(map(lambda x: int(x), line[3:])))
                # do not consider positions with coverage <= minimum
                if coverage > minimum_coverage:
                    # initialize set if key not in dict
                    if position not in sample:
                        sample[position] = set()
                    # populate dict
                    sample[position].add(participant)
        # close file after reading
        infile.close()
        
    return sample
            
# use this function to get all positions with coverage
def PositionsWithCoverage(basecall, minimum_coverage):
    '''
    (file, int) -> set
    Take the file with the basecalls and return a set of positions with
    minimum coverage (ie. positions with mapped reads and read depth > minimum)
    '''
    
    # open file for reading
    infile = open(basecall, 'r')
    # create set to store the positions
    coverage = set()
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # grab position in mito 0-based
            position = int(line[1]) -1
            # compute coverage
            reads = sum(list(map(lambda x: int(x), line[3:])))
            assert type(reads) == int, "reads should be an integer"
            # do not consider positions with read depth < minimum
            if reads > minimum_coverage:
                coverage.add(position)
    infile.close()
    return coverage
            



# use this function to get informations for positions with heteroplasmies 
def GetIndividualTumorHeteroplasmies(heteroplasmy_file, sample_size, mito_annotation):
    '''
    (file, dict, dict) -> dict
    Return a dictionary with information about the positions with heteroplasmies
    and adding the sample size for each position (ie. number of individuals
    with coverage at that position, with variant or not) and the gene names in
    mito_annotation
    '''
    # create a dict with position as key and a list of items as value
    # {position: [participant_ID, reference, gene, gene_orientation, exonic_function, AA_change, ref, major, minor, major_count, minor_count, sample_size,
    # Forward_A, Forward_T, Forward_C, Forward_G, Reverse_A, REverse_T, Reverse_C, Reverse_G, Fisher_P]
    
    variants = {}    
    
    # open file for reading
    infile = open(heteroplasmy_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # populate dict
            # make position 0-based
            position = int(line[1]) -1
            # get ID from filename   
            if 'GRCH37' in heteroplasmy_file:
                participant = heteroplasmy_file[:heteroplasmy_file.index('_')]
            elif 'RNASEQ' in heteroplasmy_file:
                participant = heteroplasmy_file[:heteroplasmy_file.index('_')]
            # get gene name and orientation
            gene = line[18]
            # check gene name
            PosInGene = False
            for MTgene in mito_annotation:
                if position in mito_annotation[MTgene]:
                    # exit loop and compare gene name
                    PosInGene = True
                    break
            # compare gene names
            if gene != '':
                if gene != 'ATP8|ATP6' and gene != 'ND4L|ND4' and gene != 'TRNI|TRNQ' and gene != 'TRNC|TRNY':
                    # note that ATP6 and ATP8 are overlapping over a short region
                    # and genes ND4 and ND4L are also overlapping
                    #  and genes TRNI and TRNQ are also overlapping
                    assert gene == MTgene, '{0} and {1} are different genes'.format(gene, MTgene)
            else:
                # check if position in annotated gene, assign MTgenee to gene
                if PosInGene == True:
                    gene = MTgene
                else:
                    gene = 'NA'
            if '(+)' in line[19]:
                orientation = 'L_(+)'
            elif '(-)' in line[19]:
                orientation = 'H_(-)'
            else:
                orientation = 'NA'
            # get exonic function
            if line[20] == '':
                exonic_function = 'NA'
            else:
                exonic_function = line[20]
            # get AA effect
            if line[21] == '':
                AA_change = 'NA'
            else:
                AA_change = line[21]
            # get ref , major and minor alelles
            ref, major, minor = line[2], line[14], line[15]
            # do bit of QC
            assert major in ['A', 'C', 'T', 'G'], 'major allele not a valid base'
            assert minor in ['A', 'C', 'T', 'G'], 'minor allele not a valid base'
            # get major and minor allele counts
            major_count, minor_count = line[16], line[17]
            # get sample size
            N = str(len(sample_size[position]))
            # get read counts for each nucleotide
            forward_A, forward_T, forward_C, forward_G = line[3], line[4], line[5], line[6]
            reverse_A, reverse_T, reverse_C, reverse_G = line[7], line[8], line[9], line[10]            
            # get the P_value of the Fisher test for heteroplasmy detection (ie freq != 0)
            Pval = line[26]
            
            # populate dict (use gene name from the rCRS annotation)
            variants[position] = [participant, gene, orientation, exonic_function,\
            AA_change, ref, major, minor, major_count, minor_count, N, forward_A,\
            forward_T, forward_C, forward_G, reverse_A, reverse_T, reverse_C, reverse_G, Pval]

    # close file after reading
    infile.close()

    return variants    


# use this function to get the position of all genes in the mitochondrial genome
def MitoAnnotation(MitoGeneFile):
    '''
    (file) -> dict
    Take a file with rCRS annotation of the mitochondrial genome and return a 
    dict with gene as key and set of position as value in 0-based    
    '''
        
    # create a dict {gene: set(positions)}
    gene_coord = {}
    
    # open file for reading
    infile = open(MitoGeneFile, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get start position 0-based
            start = int(line[1]) - 1
            # get end position
            end = int(line[2])
            # get gene name
            gene = line[3]
            gene_coord[gene] = set()
            # all indices in the set are positions in gene
            for i in range(start, end):
                gene_coord[gene].add(i)
    # close file
    infile.close()
    
    return gene_coord

# use this function to get the start and end position of the mitocondrial genes
def MitoGeneCoordinates(MitoGeneFile):
    '''
    (file) -> dict
    Take a file with rCRS annotation of the mitochondrial genome and return a 
    dict with gene as key and list of gene coordinates as value in in 0-based    
    '''
    
    # Get the position indices of the mito genes
    annotation = MitoAnnotation(MitoGeneFile)

    # create a dict with strand orientation
    strand = {'ATP6': '+', 'ATP8': '+', 'COX1': '+', 'COX2': '+', 'COX3': '+',
              'CYTB': '+', 'ND1': '+', 'ND2': '+', 'ND3': '+', 'ND4': '+',
              'ND4L': '+', 'ND5': '+', 'ND6': '-', 'RNR1': '+', 'RNR2': '+',
              'TRNA': '-', 'TRNC': '-', 'TRND': '+', 'TRNE': '-', 'TRNF': '+',
              'TRNG': '+', 'TRNH': '+', 'TRNI': '+', 'TRNK': '+', 'TRNL1': '+',
              'TRNL2': '+', 'TRNM': '+', 'TRNN': '-', 'TRNP': '-', 'TRNQ': '-',
              'TRNR': '+', 'TRNS1': '-', 'TRNS2': '+', 'TRNT': '+', 'TRNV': '+',
              'TRNW': '+', 'TRNY': '-'}

    # create a new dict so that mito_annotation is not modified  [start, end, orientation]
    coordinates = {}
    for gene in annotation:
        coordinates[gene] = list(annotation[gene])
        coordinates[gene].sort()
        start, end = coordinates[gene][0], coordinates[gene][-1] + 1
        coordinates[gene] = [start, end]
        coordinates[gene].append(strand[gene])
    
    return coordinates


# use this function to convert the heterosummaryfile into a dict of dicts
def GetVariablePositions(HeteroSummaryFile):
    '''
    (file) -> dict
    Parse the summary file with heteroplasmies and return a dictionary 
    with participant ID as key and dictionnaries of positions : information pairs
    '''
    # create a dict of {participantID: {position : [information]} 
    participants = {}
    # open file for reading
    infile = open(HeteroSummaryFile, 'r')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get participant ID
            ID = line[1]
            # get position 0-based
            position = int(line[0]) -1
            # populate dict
            if ID in participants:
                participants[ID][position] = line
            else:
                participants[ID] = {}
                participants[ID][position] = line
    infile.close()
    return participants
    

# use this function to keep only participants with both RNAseq and WGS data
def RemoveUniqueParticipants(Participants_RNA, Participants_WGS):
    '''
    (dict, dict) -> dict, dict
    Take the dictionnaries from parsing the summary heteroplasmy files for RNAseq
    and WGS and return modified dictionnaries that include only participants
    with both RNAseq and WGS data
    '''
    
    # create a set of participant ID common to RNA and WGS datasets
    CommonID = set(ID for ID in Participants_RNA if ID in Participants_WGS)
    
    # remove participants that are unique to RNAseq and WGS datasets
    to_delete = [ID for ID in Participants_RNA if ID not in CommonID]
    for ID in to_delete:
        del Participants_RNA[ID]
        
    to_delete = [ID for ID in Participants_WGS if ID not in CommonID]
    for ID in to_delete:
        del Participants_WGS[ID]
    
    return Participants_RNA, Participants_WGS

# use this function to remove positions with coverage < threshold
def RemovePositionWithoutCoverage(Participants_RNA, suffix, minimum_coverage):
    '''
    (dict, str, int) -> dict
    Take the dictionnary of the parsed RNA summary heteroplasmy file, the suffix
    in the subfolder name with mitoseek output and return a modified dictionary
    without the positions that do not have minimum coverage in the mitoseek output for WGS    
    Precondition: The participants without WGS data have been removed    
    '''
    
    # loop over participant in dict
    for ID in Participants_RNA:
        # get the subfolder corresponding to to participant
        subfolder = ID + suffix + '/'
        # get the set of positions with coverage in WGS mitoseek output
        valid_positions = PositionsWithCoverage(subfolder + 'mito1_basecall.txt', minimum_coverage)
        # create a list of positions without coverage to remove
        to_delete = [position for position in Participants_RNA[ID] if position not in valid_positions]
        for position in to_delete:
            del Participants_RNA[ID][position]
    
    # remove participants without any positions
    to_delete = [ID for ID in Participants_RNA if len(Participants_RNA[ID]) == 0]
    for ID in to_delete:
        del Participants_RNA[ID]
        
    return Participants_RNA
        
    
    
# use this function to compute at a given site for a given individual
def PolymorphismInformationContent(read_counts):
    '''
    (list) -> float
    Return the polymorphism information content (PIC) at a given site
    '''
    # PIC = 1 - (sum of pi**2 with i from 1 to n, n = number of alleles)
    
    # get total number of reads
    total = sum(read_counts)
    # set counter variable
    PIC = 0
    # loop over number of reads for each nucleotides
    for i in range(len(read_counts)):
        # update pic
        PIC += (read_counts[i] / total)**2
    
    return 1 - PIC
    

# use this function to compute PIC for all individuals at each positions
def GenomePositionsPic(HeteroSummaryFile):
    '''
    (file) -> dict
    Return a dictionary with position as key and a list of PIC values for each
    individual 
    '''

    # create a dict with positions in MT genome 0-based
    polymorphism = {}
    # initialize dict with position indices for rCRS genome sequence
    for i in range(16569):
        polymorphism[i] = []
    
    # open file for reading
    infile = open(HeteroSummaryFile, 'r')
    # skip header
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get position 0-based
            position = int(line[0]) - 1
            # count the number of reads for each nucnucleotide
            A = int(line[12]) + int(line[16])
            T = int(line[13]) + int(line[17])
            C = int(line[14]) + int(line[18])
            G = int(line[15]) + int(line[19])
            # make a list of counts
            reads = [A, T, C, G]
            # compute PIC
            PIC = PolymorphismInformationContent(reads)
            # add PIC value for that participant to list
            polymorphism[position].append(PIC)
    
    infile.close()
    return polymorphism
    
    

# use this function to determine strand bias for snps in a sample
def StrandBias(HeteroSummaryFile):
    '''
    (file) -> dict
    Return a dictionnary with mutation count for each possible change on each strand
    Count the number of different mutations, ie. do not count more than once the
    same mutation present in multiple individuals
    '''
    
    # create a dict to record the direction of mutation
    # {L_C>T: [position]}
    MutationBias = {}
    
    # loop over file
    infile = open(HeteroSummaryFile, 'r')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line !='':
            line = line.split('\t')
            # get position 0-based
            position = int(line[0]) - 1
            # get ref, major and minor
            ref, major, minor = line[7], line[8], line[9]
            # find the direction of the change and the strand
            ReadCounts = {'FA': int(line[12]), 'FT': int(line[13]),
                          'FC': int(line[14]), 'FG': int(line[15]),
                          'RA': int(line[16]), 'RT': int(line[17]),
                          'RC': int(line[18]), 'RG': int(line[19])}
            # create a list with counts for the minor allele
            MinorCounts = []
            for base in ReadCounts:
                if base[-1] == minor:
                    MinorCounts.append([ReadCounts[base], base])
            # find the strand with the most reads for the minor allele
            MinorCounts.sort()
            # get the strand and direction of change
            SNP = MinorCounts[-1][1][0] + '_' + major + '>' + minor
            # populate dict
            # check if mutation already recorded at that position
            if SNP in MutationBias:
                if position not in MutationBias[SNP]:
                    MutationBias[SNP].append(position)
            else:
                MutationBias[SNP] = [position]
                
    infile.close()
    return MutationBias
                


# use this function to compute the median read depth for a single individual
def GetReadDepth(BaseCallFile):
    '''
    (file) -> dict
    Return a dictionary with position as key and the read depth at that position
    Precondition: use rCRS mito ref genome    
    '''
    
    # create a dict {position: number of reads}
    ReadDepth = {}
    # initialize dict
    for i in range(16569):
        ReadDepth[i] = 0
    
    # open basecall file
    infile = open(BaseCallFile, 'r')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # sum the read counts for all nucleotides
            reads = sum(list(map(lambda x:int(x), line[3:])))
            # get position 0-based
            position = int(line[1]) -1
            # update read count
            ReadDepth[position] += reads

    infile.close()
    return ReadDepth        



# use this function to compute the median read depth per individual
def MedianReadDepthIndividual(folders):
    '''
    (str) -> dict
    Loop over the Mitoseek output folders and compute the median read depth 
    over all positions for each individual
    Precondition: Folders have been filtered to keep only non-redundant 
    individual mapped with GRCH37
    '''

    # create a dict {participant: median read depth}
    ReadDepth = {}
    # loop over subfolders
    for i in folders:
        # get the participant ID
        participant = i[:i.index('_')]
        # get the position-read depth for that participant
        reads = GetReadDepth(i + '/mito1_basecall.txt')
        # get the read counts
        ReadCounts = [reads[j] for j in reads]
        # median read counts
        median_reads = np.median(ReadCounts)
        # populate dict
        ReadDepth[participant] = median_reads    
    return ReadDepth


# use this function to generate a set of blacklisted individuals to ignore in all analyses
def BlackListed(blacklisted_individual_file):
    '''
    (file) -> set
    Take the file with individuals blacklisted in the PCAWG release and
    return a set of individual IDs to ignore in all analyses
    '''
    
    # create a set of individual IDs
    blacklisted = set()    
    # open file for reading
    infile = open(blacklisted_individual_file)
    # loop over file
    for line in infile:
        # skip comment lines
        if not line.startswith('#'):
            line = line.rstrip()
            if line != '':
                # add ID to set
                blacklisted.add(line)
    infile.close()
    return blacklisted
    
    
# use this function to identify multiallelic snps
def IdentifyMultiAllelicPositions(HeteroplasmyFile, threshold):
    '''
    (file, num) -> dict
    Return a dictionary with individual ID as key and an inner dictionary
    of position: list of alleles with read numbers > threshold as value
    '''
    
    # create a dict of dict {individual : {position : [allele1, allele2, allele3, allele4]}}
    MultiAlleles = {}
    
    # open file, extract multi allelic positions
    infile = open(HeteroplasmyFile, 'r')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get position 0-based
            position = int(line[0]) -1
            # get individual ID
            individual = line[1]
            # get nucleotide reads
            reads = {'A': int(line[12]) + int(line[16]),
                     'T': int(line[13]) + int(line[17]),
                     'C': int(line[14]) + int(line[18]),
                     'G': int(line[15]) + int(line[19])}
            # get the base with maximum reads
            for base in reads:
                if reads[base] == max([i for i in reads.values()]):
                    MaxBase = base
            # check that Maxbase is major allele
            assert line[14].upper() == MaxBase.upper(), 'allele with highest reads is not the major allele'
            # get total reads
            total_reads = 0
            for base in reads:
                total_reads += reads[base]
            # compare reads to threshold
            for base in reads:
                # keep allele if number of reads > threshold
                if (reads[base] / total_reads) * 100 > threshold:
                    # check if individual is already recorded
                    if individual not in MultiAlleles:
                        # inititialze inner dict
                        MultiAlleles[individual] = {}
                    # check if position is already
                    if position in MultiAlleles[individual]:
                        MultiAlleles[individual][position].append(base)
                    else:
                        MultiAlleles[individual][position] = [base]
    infile.close()
    return MultiAlleles
                        
            
# use this function to convert SNP indices on genome to indices on genes
def GenomicPositionToGenePosition(snp_start, gene_start, gene_end, orientation):
    '''
    (int, int, int, str) -> int
    Return the position of a SNP in a gene given its position on chromosome
    Precondition: all coordinates are 0-index based    
    '''
    
    if orientation == '+':
        start = snp_start - gene_start
    elif orientation == '-':
        start = (gene_end - 1) - snp_start
    return start


# use this function to find the most 5' nonsense mutation
def FindFistStopCodon(HeteroSummaryFile, MitoGeneFile):
    '''
    (file, file) -> dict
    Take the Summary file with heteroplasmies detected by MitoSeek and the file with
    mictochondrial coordinates and return a dictionary with protein-coding gene
    name as key and the 5' most upstream position of a nonsense allele in that gene
    '''
    
    # get the gene coordinates {gene: [start, end, orientation]}
    coords = MitoGeneCoordinates(MitoGeneFile)
    
    # create dict {gene: 5' position}
    stops = {}

    # open file for reading, skip header
    infile = open(HeteroSummaryFile)
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # check if stop codon gain
            if line[4] == 'stopgain':
                # gene gene name and SNP position, 0-based
                gene = line[2]
                position = int(line[0]) -1
                # check if gene in dict
                if gene in stops:
                    # check if position is most upstream (lower for + and higher for -)
                    if coords[gene][-1] == '+':
                        assert '(+)' in line[3], 'gene should be on + strand'
                        if position < stops[gene]:
                            stops[gene] = position
                    elif coords[gene][-1] == '-':
                        assert '(-)' in line[3], 'gene should be on - strand'
                        if position > stops[gene]:
                            stops[gene] = position
                else:
                    stops[gene] = position
    infile.close()
    
    # get the positions relative to gene start in 5' orientation    
    for gene in stops:
        start = GenomicPositionToGenePosition(stops[gene], coords[gene][0], coords[gene][1], coords[gene][-1])
        # drop orientation, keep only position
        stops[gene] = start

    return stops


################### edit below





def is_PTC_positions_uniform(PTC_file, window):
    '''
    (file, int) -> None
    Performs a chi-square test of uniform distribution of the most 5' PTC along the CDS
    '''

    from scipy import stats
    range_counts = PTC_positions_along_CDS(PTC_file, window)
        
    # performs chi square test
    # note that the degree of freedom = k -1 - ddof
    # with k = number of observed frequencies
    # need to specify the parameter ddof to get the appropriate degree of freedom

    chisquare_test = stats.chisquare(range_counts)
    print('chisquare: %6.4f' % chisquare_test[0])
    print('p_value: {0}'.format(chisquare_test[1]))


def histogram_PTC_positions(PTC_file, window):
    '''
    (file, int) -> list
    Return a list of porportions of the 5' most PTC along the CDS in bins of size window
    to be use to make a histogram plot
    '''
    # count the number of PTC in each bin, taling into consideration only the most 5' PTC
    range_counts = PTC_positions_along_CDS(PTC_file, window)

    # count the number of such PTC
    stops = open(PTC_file, 'r')
    header = stops.readline()
    genes = set()
    for line in stops:
        line = line.rstrip()
        if line != '':
            line  = line.split()
            gene = line[0]
            genes.add(gene)
    total = len(genes)

    # calculate proportions:
    for i in range(len(range_counts)):
        range_counts[i] = range_counts[i] / total

    stops.close()
    return range_counts



def partition_CDS_length(PTC_file):
    '''
    (file) -> list
    Sort the genes in PTC_file according to their length
    and return a list of dictionnaries for each length quartile
    and the header of the file as last item
    '''

    # stote the PTC info in dictionnary
    stops = open(PTC_file, 'r')
    PTC = {}
    header = stops.readline()
    for line in stops:
        line = line.rstrip()
        if line != '':
            line = line.split()
            PTC[line[0]] = line

    # compute the length quartiles
    CDS_length = []
    for gene in PTC:
        CDS_length.append(int(PTC[gene][3]))

    quartiles = compute_quartiles(CDS_length)

    # sort genes into dictionnary of length quartile
    Q1, Q2, Q3, Q4 = {}, {}, {}, {}
    for gene in PTC:
        size = int(PTC[gene][3])
        if size < quartiles[0]:
            Q1[gene] = PTC[gene]
        elif size >= quartiles[0] and size < quartiles[1]:
            Q2[gene] = PTC[gene]
        elif size >= quartiles[1] and size < quartiles[2]:
            Q3[gene] = PTC[gene]
        else:
            Q4[gene] = PTC[gene]

    stops.close()
    return [Q1, Q2, Q3, Q4, header]


def save_quartiles_to_file(PTC_file):
    '''
    Sort the genes according to their CDS length, store them in dictionnaries
    of quartile length and save each dictionnary to a separate file
    '''

    quartiles = partition_CDS_length(PTC_file)
    header = quartiles[-1]

    for i in range(len(quartiles)-1):
        Q_file = open('PTC_polym_SNPS_only_Q' + str(i+1) +'.txt', 'w')
        Q_file.write(header)
        for gene in quartiles[i]:
            for item in quartiles[i][gene][:-1]:
                Q_file.write(item + '\t')
            Q_file.write(quartiles[i][gene][-1] + '\n')
        Q_file.close()

