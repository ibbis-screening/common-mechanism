import pandas as pd
import numpy as np
import os
import taxoniq
import matplotlib.cm as cm

# check BLAST output files for each query to see if they exist and have any hits
def checkfile(filename):
    if os.path.isfile(filename) == False:
        print(filename, "does not exist")
        return 0
    file = open(filename).read()
    lines = 0
    for line in file.splitlines():
        if line.startswith('#') == False:
            lines = lines + 1
#    print("\n\nLines counted: " + str(lines) + "\n\n")
    if lines == 0:
#        print(filename, "has no hits")
        return 2        
    if len(file) == 0:
        print(filename, "is empty")
        return 0
    else:
        return 1

# assign a nice colour scale to numeric values associated with taxon bins - one colour for regulated orgs, another for non-regulated
def colourscale(reg_status, counts, averages):
    rmap = cm.get_cmap('OrRd', 100)
    nrmap = cm.get_cmap('Blues', 100)
    colours = []
    for i in range(len(reg_status)):
        colours.append('rgb' + str((np.array(nrmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')')) if reg_status[i] != counts[i] else colours.append('rgb' + str((np.array(rmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')'))
    return colours


def taxdist(blast, reg_ids, query):
    # create a new row for each taxon id in a semicolon-separated list, then delete the original row with the concatenated taxon ids
    # blast here is a dataframe of blast results
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
    
    # checks which individual lines contain regulated pathogens
    cut = []
    for x in range(0, blast.shape[0]):
        try:
            t = taxoniq.Taxon(blast['subject tax ids'][x])
            # taxoniq starts ranked_lineage at the species or genus level, so check the strain taxID first
            if int(blast['subject tax ids'][x]) in set(reg_ids[0]):
                blast.loc[x,'regulated'] = True
            for level in t.ranked_lineage:
                if int(level.tax_id) in set(reg_ids[0]):
                    blast.loc[x,'regulated'] = True
                if blast.columns.str.contains(level.rank.name).any():
                    blast.loc[x,level.rank.name] = level.scientific_name
                else:
                    blast[level.rank.name] = ""
                    blast.loc[x,level.rank.name] = level.scientific_name
        except:
            print('Taxon id', blast['subject tax ids'][x], 'is missing from taxoniq database')
            cut.append(x)
#
    blast.drop(cut)
    
    # identifies regulated status of each gene
#    for gene in set(blast['subject acc.']):
#        regulated = 0
#        nonregulated = 0
#        if (blast['regulated'][blast['subject acc.'] == gene] == True).any():
#            print("Reg found")
#            regulated = 1
#        if (blast['regulated'][blast['subject acc.'] == gene] == False).any():
#            print("Non-reg found")
#            nonregulated = 1
#        if (regulated == 1 & nonregulated == 0):
#            blast['regulated'][blast['subject acc.'] == gene] = "Regulated"
#        if (regulated == 1 & nonregulated == 1):
#            blast['regulated'][blast['subject acc.'] == gene] = "Mixed"
#        if (regulated == 0 & nonregulated == 0):
#            blast['regulated'][blast['subject acc.'] == gene] = "Non-regulated"
    
    blast = blast.sort_values(by=['% identity'], ascending=False)
    
    # simplify output by putting all single instances of a species in an 'other' category
    singletons = blast.species.value_counts().index[blast.species.value_counts()==1]
    blast['species_simplified'] = blast['species']
    blast.loc[blast.species.isin(singletons), 'species_simplified'] = 'other'
    
    blast = blast.reset_index(drop=True)
#    blast.to_csv(query + ".reg_path_coords.csv", columns=["q. start", "q. end"], header=None)
    
    return blast


# read in HMMER output files
def readhmmer(fileh):
    # read through the file with comments in to find column headers
    # some column names are separated by spaces, don't know how to get around this
#     file = open(fileh, "r")
#     columns = []
#     for line in file:
#         if 'target name' in line:
#             lines = line.split("  ")
#             columns = list(filter(None, lines))
#         if columns == []:
#             print("ERROR: Failed to parse column IDs")
    columns = ['target name', 'accession','tlen', 'query name',' accession','qlen','E-value','score','bias','hit #','of','c-Evalue','i-Evalue','score2','bias','hmm from','hmm to','ali from','ali to','env from','env to','acc', 'description of target']

    hmmer = pd.read_csv(fileh, sep="\s+", comment='#', header=None, engine='python', index_col=None)
    hmmer.columns = columns
    return hmmer

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
    

# read in BLAST files and pre-format the data frame with essential info
def readblast(fileh):
    blast = pd.read_csv(fileh, sep='\t', comment='#', header=None)
    
    # read through the file with comments in to find column headers
    file = open(fileh, "r")
    columns = []
    # query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score      in real title
    # qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send      in command
    for line in file:
        if 'Fields' in line:
            fields = line.replace('# Fields: ', '').rstrip()
            columns = fields.split(', ')
    if columns == []:
        print("ERROR: Failed to parse column IDs")
    
    blast.columns = columns
    blast = blast.sort_values(by=['% identity'])
    blast['log evalue'] = -np.log10(pd.to_numeric(blast['evalue'])+1e-300)
    blast['q. coverage'] = abs(blast['q. end']-blast['q. start'])/blast['query length'].max()
    blast['s. coverage'] = abs(blast['s. end']-blast['s. start'])/blast['subject length']
    blast = blast[blast['subject tax ids'].notna()]
    blast = blast.reset_index(drop=True)
    return blast

# read in DIAMOND files and pre-format the data frame with essential info
def readdmnd(fileh):
    diamond = pd.read_csv(fileh, sep='\t', comment='#', header=None)
    columns = ['query acc.', 'subject title', 'subject acc.', 'subject tax ids', 'evalue', 'bit score', '% identity', 'query length', 'q. start', 'q. end', 'subject length', 's. start', 's. end']
    
    diamond.columns = columns
    diamond = diamond.sort_values(by=['% identity'], ascending=False)
    diamond['log evalue'] = -np.log10(pd.to_numeric(diamond['evalue'])+1e-300)
    diamond['q. coverage'] = abs(diamond['q. end']-diamond['q. start'])/diamond['query length'].max()
    diamond['s. coverage'] = abs(diamond['s. end']-diamond['s. start'])/diamond['subject length']
    
    diamond = diamond[diamond['subject tax ids'].notna()]
    diamond = diamond.reset_index(drop=True)
    
    return diamond

# trim BLAST results to the most interesting ones
def trimblast(blast):
    # rank hits by PID, if any multispecies hits contain regulated pathogens, put the regulated up top
    if 'regulated' in blast:
        blast = blast.sort_values(by=['regulated'], ascending=False)
#    print(blast[['query acc.', 'subject title', 'subject tax ids', 'regulated']])
    blast = blast.drop_duplicates(subset=['query acc.', 'q. start', 'q. end'], keep='first', ignore_index=True)
    blast = blast.sort_values(by=['% identity', 'bit score'], ascending=False)
    
    drop = []
    blast2 = blast
    # only keep  top ranked hits that don't overlap
    for query in blast['query acc.'].unique():
        df = blast[blast['query acc.'] == query]
#        print("DF for trimming:\n\n")
#        print(df)
        for i in df.index:
            for j in df.index[(i+1):]:
                # if the beginning and end of the higher rank hit both extend further than the beginning and end of the lower ranked hit, discard the lower ranked hit
                if (df.loc[i,'q. start'] <= df.loc[j,'q. start'] and df.loc[i,'q. end'] >= df.loc[j,'q. end']):
                    if j in blast2.index:
                        if df.loc[i,'subject tax ids']!="32630": # don't trim a lower ranked hit if the higher ranked hit is a synthetic construct
                            blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)
    return blast2

