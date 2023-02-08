#! /usr/bin/env python

##############################################################################
#utils.py provides functions used across multiple python scripts in the 
#CommonMechanism package 
#
#Copyright (C) 2022-2023 NTI|Bio
#This file is part of the CommonMechanism
##############################################################################
#Included functions:
#   check_blastfile(filename): checks BLAST output files for each query to check if 
#       results exist and if hits exist
#   
##############################################################################

import pandas as pd
import numpy as np
import os
import taxoniq
import matplotlib.cm as cm
import re
import sys

##############################################################################
#check_blastfile 
#usage: check BLAST output files for each query to see if they exist and have any hits
#input: 
#   - name of BLAST output file  
def check_blastfile(filename):
    if os.path.isfile(filename) == False:
        # sys.stderr.write("\t...%s does not exist\n" % filename)
        return 0
    curr_file = open(filename,'r')
    num_hits = 0 
    num_lines = 0
    for line in curr_file:
        if (len(line) > 0) and (line[0] != "#"):
            num_hits += 1
        num_lines += 1
    #Return based on file results
    if num_lines == 0:
        # sys.stderr.write("\t...%s is empty\n" % filename)
        return 0
    elif num_hits == 0:
        # sys.stderr.write("\t...%s has no hits\n" % filename)
        return 2        
    else:
        return 1

##############################################################################
#colourscale 
#usage: assigns colour scale to numeric values associated with taxon bins
#   separate colours for regulated orgs (OrRd) and non-regulated (Blues)
#input:
#   - reg_status: regulatory status
#   - counts: counts to evaluate
#   - averages: averages to evaluate 
def colourscale(reg_status, counts, averages):
    rmap = cm.get_cmap('OrRd', 100)
    nrmap = cm.get_cmap('Blues', 100)
    colours = []
    for i in range(len(reg_status)):
        colours.append('rgb' + str((np.array(nrmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')')) if reg_status[i] != counts[i] else colours.append('rgb' + str((np.array(rmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')'))
    return colours

##############################################################################
#split_taxa 
def split_taxa(blast):
    blast2 = blast
    cutrows = []
    lastrow = blast.shape[0]
    for i in range(0, len(blast['subject tax ids'])):
        if str(blast.loc[i,'subject tax ids']).find(";") != -1:
            taxids = str(blast.loc[i,'subject tax ids']).split(";")
            cutrows.append(i)
            for tax in taxids:
                blast2.loc[lastrow+1,:] = blast.loc[i,:]
                blast2.loc[lastrow+1,'subject tax ids'] = tax
                lastrow = lastrow+1
    
    blast = blast2.drop(cutrows)
    blast = blast.reset_index(drop=True)
    blast['regulated'] = False

    return blast

##############################################################################
#taxdist 
def taxdist(blast, reg_ids, vax_ids):
    # create a new row for each taxon id in a semicolon-separated list, then delete the original row with the concatenated taxon ids
    # blast here is a dataframe of blast results
    blast = split_taxa(blast)
    blast = blast[blast['subject tax ids'] != "32630"]
    blast = blast.reset_index(drop=True)
    blast['subject tax ids'] = blast['subject tax ids'].astype('int')
    # print(blast)
    
    # checks which individual lines contain regulated pathogens
    missing = []
    vax = []
    columns = []

    for x in range(0, blast.shape[0]):
        try:
            t = taxoniq.Taxon(int(blast['subject tax ids'][x]))
            # taxoniq starts ranked_lineage at the species or genus level, so check the strain taxID first
            if int(blast['subject tax ids'][x]) in set(reg_ids[0]):
                blast.loc[x,'regulated'] = True
            while t.scientific_name != 'root': #  & int(t.tax_id) not in set(vax_ids[0])
                if int(t.tax_id) in set(reg_ids[0]):
                    blast.loc[x,'regulated'] = True
                if blast.columns.str.contains(t.rank.name).any():
                    blast.loc[x,t.rank.name] = t.scientific_name
                else:
                    blast[t.rank.name] = ""
                    blast.loc[x,t.rank.name] = t.scientific_name
                    columns.append(t.rank.name)
                t = t.parent
                if int(t.tax_id) in set(vax_ids[0]):
                    vax.append(x)
                # print(t.rank.name)
                # print(t.scientific_name)
        except:
            sys.stderr.write('\t...taxon id ' + str(blast['subject tax ids'][x]) + ' is missing from taxoniq database\n')
            missing.append(x)
    
    if (len(vax)>0):
        blast['regulated'][vax] = False
#
    blast = blast.drop(missing)
    # if (len(missing)>0):
    #     blast.loc[missing,columns] = "missing"
        
    blast = blast.sort_values(by=['% identity'], ascending=False)
    # print(blast)
    
    # simplify output by putting all single instances of a species in an 'other' category
    singletons = blast.species.value_counts().index[blast.species.value_counts()==1]
    blast['species_simplified'] = blast['species']
    blast.loc[blast.species.isin(singletons), 'species_simplified'] = 'other'
    
    blast = blast.reset_index(drop=True)
    
    return blast


#####################################################################################
#readhmmer
#usage:  read in HMMER output files
#input:
#   - hmmer output file name
def readhmmer(fileh):
    columns = ['target name', 'accession','tlen', 'query name',' accession','qlen','E-value','score','bias','hit #','of','c-Evalue','i-Evalue','score2','bias','hmm from','hmm to','ali from','ali to','env from','env to','acc', 'description of target']
    
    hmmer = []

    with open(fileh, 'r') as f:
        for line in f:
            if "#" in line:
                continue
            bits = re.split('\s+', line)
            description = " ".join(bits[22:])
            bits = bits[:22]
            bits.append(description)
            hmmer.append(bits)
    hmmer = pd.DataFrame(hmmer, columns=columns)
    hmmer['E-value'] = pd.to_numeric(hmmer['E-value'])
    hmmer['score'] = pd.to_numeric(hmmer['score'])
    hmmer['ali from'] = pd.to_numeric(hmmer['ali from'])
    hmmer['ali to'] = pd.to_numeric(hmmer['ali to'])
    hmmer['qlen'] = pd.to_numeric(hmmer['qlen'])

#    print(hmmer)
    return hmmer

#####################################################################################
def trimhmmer(hmmer): # don't forget this is a report on 6-frame translations so coordinates will be complicated
    # rank hits by bitscore
    hmmer = hmmer.sort_values(by=['score'], ascending=False)
#     hmmer = hmmer.drop_duplicates(subset=['query acc.', 'q. start', 'q. end'], keep='first', ignore_index=True)
    
    drop = []
    hmmer2 = hmmer
    # only keep  top ranked hits that don't overlap
    for query in hmmer['query name'].unique():
        df = hmmer[hmmer['query name'] == query]
        for i in df.index:
            for j in df.index[(i+1):]:
                if (df.loc[i,'ali from'] <= df.loc[j,'ali from'] and df.loc[i,'ali to'] >= df.loc[j,'ali to']) | (df.loc[i,'ali from'] >= df.loc[j,'ali from'] and df.loc[i,'ali to'] <= df.loc[j,'ali to']):
                    if j in hmmer2.index:
                            hmmer2 = hmmer2.drop([j])
        hmmer2 = hmmer2.reset_index(drop=True)
    return hmmer2
    

#####################################################################################
# read in BLAST/DIAMOND files and pre-format the data frame with essential info
def readblast(fileh):
    blast = pd.read_csv(fileh, sep='\t', comment='#', header=None)
    columns = ['query acc.', 'subject title', 'subject acc.', 'subject tax ids', 'evalue', 'bit score', '% identity', 'query length', 'q. start', 'q. end', 'subject length', 's. start', 's. end']
    
    blast.columns = columns
    blast = blast.sort_values(by=['% identity'], ascending=False)
    blast['log evalue'] = -np.log10(pd.to_numeric(blast['evalue'])+1e-300)
    blast['q. coverage'] = abs(blast['q. end']-blast['q. start'])/blast['query length'].max()
    blast['s. coverage'] = abs(blast['s. end']-blast['s. start'])/blast['subject length']
    
    blast = blast[blast['subject tax ids'].notna()]
    blast = blast.reset_index(drop=True)
    
    return blast

#####################################################################################
# trim BLAST results to the most interesting ones
def trimblast(blast):
    # set start to the lowest coordinate and end to the highest to avoid confusion
    for j in blast.index:
        if blast.loc[j,'q. start'] > blast.loc[j,'q. end']:
            start = blast.loc[j,'q. end']
            end = blast.loc[j,'q. start']
            blast.loc[j,'q. start'] = start
            blast.loc[j,'q. end'] = end
            # print(blast.loc[j,['subject acc.', 'q. start', 'q. end', 's. start', 's. end']])

    # rank hits by PID, if any multispecies hits contain regulated pathogens, put the regulated up top
    if 'regulated' in blast:
        blast = blast.sort_values(by=['regulated'], ascending=False)
    blast = blast.drop_duplicates(subset=['query acc.', 'q. start', 'q. end'], keep='first', ignore_index=True)
    blast = blast.sort_values(by=['% identity', 'bit score'], ascending=False)
    blast = blast.reset_index(drop=True)
    
    blast2 = blast
    # only keep  top ranked hits that don't overlap
    for query in blast['query acc.'].unique():
        df = blast[blast['query acc.'] == query]
        # print(df)
        for i in df.index: # run through each hit from the top
            for j in df.index[(i+1):]: # compare to each below
                # if the beginning and end of the higher rank hit both extend further than the beginning and end of the lower ranked hit, discard the lower ranked hit
                if (df.loc[i,'q. start'] <= df.loc[j,'q. start'] and df.loc[i,'q. end'] >= df.loc[j,'q. end']):
                    if j in blast2.index:
                        if df.loc[i,'subject tax ids']!="32630": # don't trim a lower ranked hit if the higher ranked hit is a synthetic construct
                            if (df.loc[i,'q. start'] < df.loc[j,'q. start'] or df.loc[i,'q. end'] > df.loc[j,'q. end'] or df.loc[i,'% identity'] > df.loc[j,'% identity']): # don't drop hits if they have the same coordinates and % identity
                                blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)
#    print(blast2[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end']])
    
    return blast2

# go through trimmed BLAST hits and only look at top protein hit for each base
def tophits(blast2):
    print(blast2.loc[:10,['subject acc.', 'q. start', 'q. end', 's. start', 's. end']])
    
    # set start to the lowest coordinate and end to the highest to avoid confusion
    # for j in blast2.index:
    #     if blast2.loc[j,'q. start'] > blast2.loc[j,'q. end']:
    #         start = blast2.loc[j,'q. end']
    #         end = blast2.loc[j,'q. start']
    #         blast2.loc[j,'q. start'] = start
    #         blast2.loc[j,'q. end'] = end
    #         print(blast2.loc[j,['subject acc.', 'q. start', 'q. end', 's. start', 's. end']])
    
    blast3 = blast2
    # only keep coordinates of each hit that are not already covered by a better hit
    for query in blast3['query acc.'].unique():
        df = blast3[blast3['query acc.'] == query]
        # print(df)
        for i in df.index: # run through each hit from the top
            for j in df.index[(i+1):]: # compare to each below
                # if the beginning of a weaker hit is inside a stronger hit, alter its start to the next base after that hit
                if (df.loc[j,'q. start'] >= df.loc[i,'q. start'] and df.loc[j,'q. start'] <= df.loc[i,'q. end']):
                    # print(df.loc[j,'subject acc.'], df.loc[j,'q. start'], df.loc[j,'q. end'], df.loc[i,'q. start'], df.loc[i,'q. end'], df.loc[i,'q. end'] + 1)
                    if (df.loc[i,'q. end'] + 1 < df.loc[j,'q. end']):
                        df.loc[j,'q. start'] = df.loc[i,'q. end'] + 1
                    else:
                        df.loc[j,'q. start'] = 0
                        df.loc[j,'q. end'] = 0
                # if the end of a weaker hit is inside a stronger hit, alter the end to just before that hit
                if (df.loc[j,'q. end'] >= df.loc[i,'q. start'] and df.loc[j,'q. end'] <= df.loc[i,'q. end']):
                    # print(df.loc[j,'subject acc.'], df.loc[j,'q. end'], df.loc[i,'q. start'], df.loc[i,'q. end'], df.loc[i,'q. start'] - 1)
                    if (df.loc[i,'q. end'] + 1 < df.loc[j,'q. start']):
                        df.loc[j,'q. end'] = df.loc[i,'q. start'] - 1
                    else:
                        df.loc[j,'q. start'] = 0
                        df.loc[j,'q. end'] = 0
        for j in df.index: 
            blast3.loc[j,'subject length'] = max([df.loc[j,'q. start'], df.loc[j,'q. end']]) - min([df.loc[j,'q. start'], df.loc[j,'q. end']])
            blast3.loc[j,'q. start'] = df.loc[j,'q. start']
            blast3.loc[j,'q. end'] = df.loc[j,'q. end']
                
    blast3 = blast3[blast3['subject length'] >= 50]
    # blast3 = blast3[blast3['q. start'] < blast3['q. end']]
    print(blast3[['subject acc.', 'q. start', 'q. end', 's. start', 's. end']][:40])
    print(blast2.loc[blast3.index,['subject acc.', 'q. start', 'q. end', 's. start', 's. end']][:40])
    blast3 = blast3.reset_index(drop=True)
    # print(blast3['subject length'])
    return blast3