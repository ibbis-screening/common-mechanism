from utils import *
import os, sys
import pandas as pd
import taxoniq

# check for any best matches to a regulated pathogen in the BLAST results, and if so, print their coordinates

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

file = sys.argv[1] + ".nr.blastx"
reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)

# read in BLAST output and label with regulated pathogens
def taxdist(file, reg_ids):
    if checkfile(file) == 0:
        return
    blast = readblast(file)
    
    # create a new row for each taxon id in a semicolon-separated list, then delete the original row with the concatenated taxon ids
    blast2 = blast
    cutrows = []
    lastrow = len(blast['subject tax ids'])
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
    
    # for each row, check taxon IDs to see if they fall under any regulated pathogen taxons
    for x in range(0, blast.shape[0]):
        try:
            t = taxoniq.Taxon(blast['subject tax ids'][x])
            for level in t.ranked_lineage:
                if int(level.tax_id) in set(reg_ids[0]): # changed reg_ids to set(reg_ids[0])
                    blast.loc[x,'regulated'] = True
                if blast.columns.str.contains(level.rank.name).any():
                    blast.loc[x,level.rank.name] = level.scientific_name
                else:
                    blast[level.rank.name] = ""
                    blast.loc[x,level.rank.name] = level.scientific_name
        except:
            print('FYI: Taxon id', blast['subject tax ids'][x], 'is missing from taxoniq records (if this is a concern, contact taxoniq developers)')
        
    blast = blast.reset_index(drop=True)
    
    return blast

blast = taxdist(file, reg_ids)
blast = trimblast(blast)

if blast['regulated'].sum():
    print("Regulated pathogens: FLAG")
    hits = blast[blast['regulated']==True][['q. start', 'q. end']]   # print out the start and end coordinated on the query sequence
    hits.to_csv(sys.argv[1] + ".reg_path_coords.csv")
else:
	print("Regulated pathogens: PASS")
