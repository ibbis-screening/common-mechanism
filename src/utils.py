import pandas as pd
import numpy as np
import os

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
    if lines == 0:
        # print(filename, "has no hits")
        return 2        
    if len(file) == 0:
        print(filename, "is empty")
        return 0
    else:
        return 1

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
    columns = ['target name', 'accession','tlen', 'query name',' accession','qlen',' E-value','score','bias','hit #','of','c-Evalue','i-Evalue','score','bias','hmm from','hmm to','ali from','ali to','env from','env to','acc', 'description of target']

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
                if (df.loc[i,'q. start'] <= df.loc[j,'q. start'] and df.loc[i,'q. end'] >= df.loc[j,'q. end']) | (df.loc[i,'q. start'] >= df.loc[j,'q. start'] and df.loc[i,'q. end'] <= df.loc[j,'q. end']):
                    if j in blast2.index:
                        blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)
    return blast2
    

# read in BLAST files and pre-format the data frame with essential info
def readblast(fileh):
    blast = pd.read_csv(fileh, sep='\t', comment='#', header=None)
    
    # read through the file with comments in to find column headers
    file = open(fileh, "r")
    columns = []
    for line in file:
        if 'Fields' in line:
            fields = line.replace('# Fields: ', '').rstrip()
            columns = fields.split(', ')
    if columns == []:
        print("ERROR: Failed to parse column IDs")
    
    blast.columns = columns
    blast.sort_values(by=['% identity'])
    blast['log evalue'] = -np.log10(pd.to_numeric(blast['evalue'])+1e-300)
    blast['q. coverage'] = abs(blast['q. end']-blast['q. start'])/blast['query length'].max()
    blast['s. coverage'] = abs(blast['s. end']-blast['s. start'])/blast['subject length']
    blast = blast[blast['subject tax ids'].notna()]
    blast = blast.reset_index(drop=True)
    return blast

# trim BLAST results to the most interesting ones
def trimblast(blast):
    # rank hits by PID
    blast = blast.sort_values(by=['% identity', 'bit score'], ascending=False)
    blast = blast.drop_duplicates(subset=['query acc.', 'q. start', 'q. end'], keep='first', ignore_index=True)
    
    drop = []
    blast2 = blast
    # only keep  top ranked hits that don't overlap
    for query in blast['query acc.'].unique():
        df = blast[blast['query acc.'] == query]
        for i in df.index:
            for j in df.index[(i+1):]:
                if (df.loc[i,'q. start'] <= df.loc[j,'q. start'] and df.loc[i,'q. end'] >= df.loc[j,'q. end']) | (df.loc[i,'q. start'] >= df.loc[j,'q. start'] and df.loc[i,'q. end'] <= df.loc[j,'q. end']):
                    if j in blast2.index:
                        blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)
    return blast2

