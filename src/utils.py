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
# import taxoniq
import matplotlib.cm as cm
import re
import sys
from tqdm import tqdm
import pytaxonkit


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
#   - counts: count of hits to evaluate
#   - averages: averages to evaluate 
def colourscale(reg_status, counts, averages):
    rmap = cm.get_cmap('OrRd', 100)
    nrmap = cm.get_cmap('Blues', 100)
    colours = []
    for i in range(len(reg_status)):
        colours.append('rgb' + str((np.array(nrmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')')) if reg_status[i] != counts[i] else colours.append('rgb' + str((np.array(rmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')'))
    return colours

##############################################################################
#split_taxa: splits multi-taxon IDs in BLAST results into multiple rows in the results table each with their own taxon ID
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
    blast['superkingdom'] = ""
    blast['species'] = ""
    blast['phylum'] = ""

    return blast

##############################################################################
#taxdist 
 
def taxdist(blast, reg_ids, vax_ids, db_path, threads):
    # prevent truncation of taxonomy results
    pd.set_option('display.max_colwidth', None)

    # create a new row for each taxon id in a semicolon-separated list, then delete the original row with the concatenated taxon ids
    # blast here is a dataframe of blast results
    blast = split_taxa(blast)
    blast['subject tax ids'] = blast['subject tax ids'].astype('int')
    blast = blast[blast['subject tax ids'] != 32630] # synthetic constructs
    blast = blast[blast['subject tax ids'] != 29278] # vectors
    blast = blast.reset_index(drop=True)
    # print(blast)
    
    # checks which individual lines contain regulated pathogens
    t = pytaxonkit.lineage(blast['subject tax ids'], data_dir=db_path, threads=threads)
    reg = list(map(str, reg_ids[0]))
    vax = list(map(str, vax_ids[0]))

    for x in range(0, blast.shape[0]): # for each hit taxID
        # fetch the full lineage for that taxID
        

        # go through each taxonomy level and check for regulated taxIDs
        tax_lin = pd.DataFrame(list(zip(t['FullLineage'].str.split(';')[x], t['FullLineageTaxIDs'].str.split(';')[x], t['FullLineageRanks'].str.split(';')[x])), columns=['Lineage', 'TaxID', 'Rank'])
        tax_lin.set_index('Rank', inplace=True)
        # print(tax_lin)

        taxlist = list(map(str, tax_lin['TaxID']))
        exlist = ['32630', '29278']    

        if any(x in exlist for x in taxlist):
            blast.drop(x, axis=0, inplace=True)
            print("Drop")
            continue
        if any(x in reg for x in taxlist):
            # blast.loc[x,'regulated'] = ",".join([x for x in taxlist if x in reg])
            blast.loc[x,'regulated'] = True
            # print("Flag")
        if any(x in vax for x in taxlist):
            blast.loc[x,'regulated'] = False

        if 'superkingdom' in tax_lin.index:
            blast.loc[x,'superkingdom'] = tax_lin.loc['superkingdom', 'Lineage']
        else:
            blast.loc[x,'superkingdom'] = ""
        if 'species' in tax_lin.index:
            blast.loc[x,'species'] = tax_lin.loc['species', 'Lineage']
        else:
            blast.loc[x,'species'] = ""
        if 'phylum' in tax_lin.index:
            blast.loc[x,'phylum'] = tax_lin.loc['phylum', 'Lineage']
        else:
            blast.loc[x,'phylum'] = ""
      
    blast = blast.sort_values(by=['% identity'], ascending=False)
  
    blast = blast.reset_index(drop=True)
  
    return blast

def taxdist_old(blast, reg_ids, vax_ids):
    # create a new row for each taxon id in a semicolon-separated list, then delete the original row with the concatenated taxon ids
    # blast here is a dataframe of blast results
    blast = split_taxa(blast)
    blast['subject tax ids'] = blast['subject tax ids'].astype('int')
    blast = blast[blast['subject tax ids'] != 32630]
    blast = blast.reset_index(drop=True)
    # print(blast)
    
    # checks which individual lines contain regulated pathogens
    missing = []
    vax = []
    artif = []
    columns = []

    for x in range(0, blast.shape[0]):
        try:
            t = taxoniq.Taxon(int(blast['subject tax ids'][x]))
            # taxoniq starts ranked_lineage at the species or genus level, so check the strain taxID first
            if int(blast['subject tax ids'][x]) in set(reg_ids[0]):
                # blast.loc[x,'regulated'] = True
                blast.loc[x,'regulated'] = int(blast['subject tax ids'][x])
            while t.scientific_name != 'root': #  & int(t.tax_id) not in set(vax_ids[0])
                if int(t.tax_id) in set(reg_ids[0]):
                    # blast.loc[x,'regulated'] = True
                    blast.loc[x,'regulated'] = int(t.tax_id)
                if int(t.tax_id) == 81077:
                    artif.append(x)
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
    blast = blast.drop(artif)
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
            if "# Program:         hmmscan" in line:
                break
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
#readcmcan
#usage:  read in cmscan output files
#input:
#   - hmmer output file name
def readcmscan(fileh):
    columns = ['target name', 'accession','query name','accession','mdl','mdl from','mdl to', 'seq from', 'seq to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description of target']
    
    cmscan = []

    with open(fileh, 'r') as f:
        for line in f:
            if "# Program:         cmscan" in line:
                break
            if "#" in line:
                continue
            bits = re.split('\s+', line)
            description = " ".join(bits[17:])
            bits = bits[:17]
            bits.append(description)
            cmscan.append(bits)
    cmscan = pd.DataFrame(cmscan, columns=columns)
    cmscan['E-value'] = pd.to_numeric(cmscan['E-value'])
    cmscan['score'] = pd.to_numeric(cmscan['score'])
    cmscan['seq from'] = pd.to_numeric(cmscan['seq from'])
    cmscan['seq to'] = pd.to_numeric(cmscan['seq to'])

#    print(cmscan)
    return cmscan

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
    blast = shift_hits_pos_strand(blast)

    # rank hits by PID, if any multispecies hits contain regulated pathogens, put the regulated up top
    if 'regulated' in blast:
        blast = blast.sort_values(by=['regulated'], ascending=False)
    blast = blast.sort_values(by=['% identity', 'bit score'], ascending=False)
    # blast = blast.drop_duplicates(subset=['query acc.', 'subject acc.', 'q. start', 'q. end'], keep='first', ignore_index=True)
    blast = blast.reset_index(drop=True)
    # print(blast[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']])

    blast2 = blast
    # only keep  top ranked hits that don't overlap
    for query in blast['query acc.'].unique():
        df = blast[blast['query acc.'] == query]
        # print(df)
        for i in df.index: # run through each hit from the top
            for j in df.index[(i+1):]: # compare to each below
                if j in blast2.index:
                    # if the beginning and end of the higher rank hit both overlap or extend further than the beginning and end of the lower ranked hit, discard the lower ranked hit
                    if (df.loc[i,'q. start'] <= df.loc[j,'q. start'] and df.loc[i,'q. end'] >= df.loc[j,'q. end']):
                        if (df.loc[i,'q. start'] < df.loc[j,'q. start'] or df.loc[i,'q. end'] > df.loc[j,'q. end'] or df.loc[i,'% identity'] > df.loc[j,'% identity']): # don't drop hits if they have the same coordinates and % identity
                            # print(str(i) + " " + str(j))
                            blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)
    # print(blast2[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']])
    
    return blast2

def trim_to_top(df):
    keep_rows = []
    df = df.sort_values('% identity', ascending=False)
    df.reset_index(inplace=True)
    df = shift_hits_pos_strand(df)
    prev_hit = None
    for base in range(1, df['query length'][0]):
        # identify the row index of the top scoring hit
        if df[(df['q. start'] <= base) & (df['q. end'] >= base).all()].shape[0] > 0:
            top_hit = df[(df['q. start'] <= base) & (df['q. end'] >= base).all()].index[0]
            # if this hit hasn't been top before, set the start of the query coverage to this base
            if top_hit not in keep_rows:
                df.loc[top_hit,'q. start'] = base
                keep_rows.append(top_hit)
            # if the top hit just changed, set the end of query coverage for the last hit to the previous base
            if (top_hit != prev_hit) & (prev_hit != None):
                df.loc[prev_hit, 'q. end'] = base - 1
            prev_hit = top_hit
    
    # print(keep_rows)
    return df.iloc[keep_rows]

def shift_hits_pos_strand(blast):
    for j in blast.index:
        if blast.loc[j,'q. start'] > blast.loc[j,'q. end']:
            start = blast.loc[j,'q. end']
            end = blast.loc[j,'q. start']
            blast.loc[j,'q. start'] = start
            blast.loc[j,'q. end'] = end
    return blast

def trim_edges(df):
    for top in range(len(df.index)): # run through each hit from the top
            i = df.index[top]
            # if [df['subject acc.'][i] == "EOB10733.1"]:
            #     print("EOB10733.1")
            for next in range(top+1, len(df.index)): # compare to each below
                j = df.index[next]
                i_start = df.loc[i,'q. start']
                i_end = df.loc[i,'q. end']
                j_start = df.loc[j,'q. start']
                j_end = df.loc[j,'q. end']
                
                # if the beginning of a weaker hit is inside a stronger hit, alter its start to the next base after that hit
                if (j_start >= i_start and j_start <= i_end):
                    # print(df.loc[j,'subject acc.'], df.loc[j,'q. start'], df.loc[j,'q. end'], df.loc[i,'q. start'], df.loc[i,'q. end'], df.loc[i,'q. end'] + 1)
                    # keep equivalent hits
                    # if (j_start == i_start and j_end == i_end and df.loc[j,'% identity'] == df.loc[i,'% identity']):
                    if (df.loc[j,'% identity'] == df.loc[i,'% identity']):
                        pass
                    # if the hit extends past the end of the earlier one
                    elif (i_end + 1 < j_end):
                        # if df['subject acc.'][j] == 'WP_235443889.1':
                        # print('Trimming ' + str(j) + " based on " + str(i))
                        df.loc[j,'q. start'] = i_end + 1
                        # print(df.loc[j,['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']])
                    elif (i_end == j_end and df.loc[j,'% identity'] == df.loc[i,'% identity']):
                        pass
                    # remove if the hit is contained in the earlier one
                    else:
                        df.loc[j,'q. start'] = 0
                        df.loc[j,'q. end'] = 0
                
                # if the end of a weaker hit is inside a stronger hit, alter the end to just before that hit
                if (j_end >= i_start and j_end <= i_end):
                    # print(df.loc[j,'subject acc.'], df.loc[i,'subject acc.'], df.loc[j,'q. end'], df.loc[i,'q. start'], df.loc[i,'q. end'], df.loc[i,'q. start'] - 1)
                    # keep equivalent hits
                    # if (j_start == i_start and j_end == i_end and df.loc[j,'% identity'] == df.loc[i,'% identity']):
                    if (df.loc[j,'% identity'] == df.loc[i,'% identity']):
                        pass
                    elif (i_start - 1 > j_start):
                        # if df['subject acc.'][j] == 'WP_235443889.1':
                        #     print('Trimming')
                        #     print(df.loc[j,'subject acc.'], df.loc[i,'subject acc.'], df.loc[j,'q. end'], df.loc[i,'q. start'], df.loc[i,'q. end'], df.loc[i,'q. start'] - 1)
                        df.loc[j,'q. end'] = i_start - 1
                    elif (i_start == j_start and df.loc[j,'% identity'] == df.loc[i,'% identity']):
                        pass
                    else:
                        df.loc[j,'q. start'] = 0
                        df.loc[j,'q. end'] = 0
    
    rerun = 0
    mix_starts = 0
    for start in set(df['q. start']):
        if len(set(zip(df['q. start'][df['q. start'] == start], df['q. end'][df['q. start'] == start]))) > 1:
            if len(set(df['% identity'][df['q. start'] == start])) > 1: # if there are overlapping annotations with different % identities, re-run
                # print("Found a length discrepancy\n\n")
                # print(df[df['q. start'] == start][['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']])
                # print(set(zip(df['q. start'][df['q. start'] == start], df['q. end'][df['q. start'] == start])))
                # print(len(set(zip(df['q. start'][df['q. start'] == start], df['q. end'][df['q. start'] == start]))))
                rerun = 1
                mix_starts = mix_starts + 1
    # if rerun == 1:
    #     print('Running again - ' + str(mix_starts) + ' non-consistent starts identified')
    return df, rerun
        

# go through trimmed BLAST hits and only look at top protein hit for each base
def tophits(blast2):
    
    blast3 = blast2
    # print(blast2)
    # blast3 = blast3.sort_values('q. start')
    blast3 = blast3.sort_values('% identity', ascending=False)

    # only keep coordinates of each hit that are not already covered by a better hit
    for query in blast3['query acc.'].unique():
        # print("Filtering to top hits")
        # print(blast3)
        df = blast3[blast3['query acc.'] == query]
        # df.reset_index(inplace=True)
        # print(df[['query acc.', 'subject title', 'subject tax ids', 'regulated', 'q. start', 'q. end', '% identity']].head(20))

        rerun = 1
        while rerun == 1: # edges of hits can be moved within a higher scoring hit in the first pass
            df, rerun = trim_edges(df)
        # df, rerun = trim_edges(df)

        for j in df.index: 
            blast3.loc[j,'subject length'] = max([df.loc[j,'q. start'], df.loc[j,'q. end']]) - min([df.loc[j,'q. start'], df.loc[j,'q. end']])
            blast3.loc[j,'q. start'] = df.loc[j,'q. start']
            blast3.loc[j,'q. end'] = df.loc[j,'q. end']
        # print("Done filtering")
        # print(blast3)
                
    # print("Results of tophit")
    blast3 = blast3.sort_values('q. start')
    blast3 = blast3[blast3['q. start'] != 0]
    # print(blast3[['subject acc.', 'q. start', 'q. end', 's. start', 's. end']][:40])
    
    # only keep annotations covering 50 bases or more
    blast3 = blast3[blast3['subject length'] >= 50]
    # blast3 = blast3[blast3['q. start'] < blast3['q. end']]
    # print(blast3[['subject acc.', 'q. start', 'q. end', 's. start', 's. end']][:40])
    # print(blast2.loc[blast3.index,['subject acc.', 'q. start', 'q. end', 's. start', 's. end']][:40])
    blast3 = blast3.reset_index(drop=True)
    # print(blast3['subject length'])
    return blast3