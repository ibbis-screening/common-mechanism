import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import plotly
import matplotlib.cm as cm
import matplotlib
import re
import os
import subprocess
from Bio import ExPASy
from Bio import SwissProt
import taxoniq

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
    
def readblast(file):
    blast = pd.read_csv(file, sep='\t', comment='#', header=None)
    file = open(file, "r")
    for line in file:
        if re. search('Fields', line):
            fields = line.replace('# Fields: ', '').rstrip()
#             colummns.append(fields.split(', '))
            columns = fields.split(', ')
    blast.columns = columns       
    blast.sort_values(by=['% identity'])
    blast['log evalue'] = -np.log10(pd.to_numeric(blast['evalue'])+1e-300)
    blast['q. coverage'] = abs(blast['q. end']-blast['q. start'])/blast['query length'].max()
    blast['s. coverage'] = abs(blast['s. end']-blast['s. start'])/blast['subject length']
    blast = blast[blast['subject tax ids'].notna()]
    blast = blast.reset_index(drop=True)
    return blast

def trimblast(blast):
    blast = blast.sort_values(by='% identity', ascending=False)
    blast = blast.drop_duplicates(subset=['q. start', 'q. end'], keep='first', ignore_index=True)
    
    drop = []
    blast2 = blast
    for i in range(blast.shape[0]):
        for j in blast.index[(i+1):]:
            if (blast.loc[i,'q. start'] <= blast.loc[j,'q. start'] and blast.loc[i,'q. end'] >= blast.loc[j,'q. end']) | (blast.loc[i,'q. start'] >= blast.loc[j,'q. start'] and blast.loc[i,'q. end'] <= blast.loc[j,'q. end']):
                if j in blast2.index:
                    blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)
    return blast2

def plothits(starts, ends, qlen, names, colours, nhits):
    if nhits <5:
        yax = 6
    else:
        yax=nhits+1
    fig = go.Figure(go.Scatter(x=[0,0], y=[0,0],mode="markers",marker=dict(color=[0,100], colorscale="YlOrRd",colorbar=dict(title="% identity", x=-0.15, xanchor="left"))), go.Layout(plot_bgcolor="white"))
    fig.add_shape(type="rect", name = 'Query', x0 = 1, x1 = qlen, y0=0.5, y1=1.3, line=dict(color="white"), fillcolor='grey')
    fig.add_annotation(xanchor='left', text='Query', x=qlen, y=(0.5+1.3)/2,font=dict(family="Arial Narrow", size=int(100/yax), color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
    for i in range(0, starts.shape[0]):
        start = starts[i]
        end = ends[i]
        hit_description = names[i]
        colour = colours[i]
        fig.add_shape(type="rect", name = hit_description, x0 = start, x1 = end, y0=0.5+i+1, y1=1.3+i+1, line=dict(color="white"), fillcolor=matplotlib.colors.rgb2hex(cm.YlOrRd(colour/100)))
        fig.add_trace(go.Scatter(text = str(colour), x = [start, end, end, start, start], y=[0.5+i+1, 0.5+i+1, 1.3+i+1, 1.3+i+1, 0.5+i+1], fill="toself", line=dict(color="white")))
        fig.add_annotation(xanchor='left', text=hit_description, x=qlen, y=(0.5+i+1+1.3+i+1)/2,font=dict(family="Arial Narrow", size=int(100/yax), color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
    fig.update_xaxes(range=[0, qlen*2+3000])
    fig.update_yaxes(range=[yax+0.5, 0.5], showticklabels=False)
    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0))
    return fig

def plot_sunburst(blast, query, suffix):
    # visualisation  
    fig = px.sunburst(blast, path=['superkingdom', 'phylum', 'class', 'genus', 'species', 'subject title'], color_continuous_scale=px.colors.sequential.YlOrRd, title="Taxonomic distribution of hits", color=-np.log10(pd.to_numeric(blast['evalue'])+1e-300), labels={"color": "-log10(E-value)"})
    fig.update_layout(margin=dict(l=20, r=20, t=30, b=0))
    fig.write_image("figures/" + query + "." + suffix + "_tax.png", width=900, height=700, scale=2)

def genSankey(df, cat_cols=[], value_cols='count', colour_col='evalue', title='Sankey Diagram'):
    # maximum of 6 value cols -> 6 colors
    colorPalette = ['#4B8BBE','#306998','#FFE873','#FFD43B','#646464', "yellow", "green", "blue"]
    xlayout = [0.1,0.2,0.3,0.4,0.7]
    labelList = []
    colorNumList = []
#     for catCol in cat_cols:
# #         print(list(set(df[catCol].values)))
#         labelListTemp =  list(set(df[catCol].values))
#         colorNumList.append(len(labelListTemp))
#         labelList = labelList + labelListTemp
        
#     # remove duplicates from labelList
#     labelList = list(dict.fromkeys(labelList))
# #     print(labelList)
    
    # define colors based on number of levels
    colorList = []
    for idx, colorNum in enumerate(colorNumList):
        colorList = colorList + [colorPalette[idx]]*colorNum
        
    # transform df into a source-target pair
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
            mean_val = [df[colour_col][df[cat_cols[i]]==x].mean() for x in df[cat_cols[i]]]
            sourceTargetDf.loc[:,'colour'] = mean_val
            sourceTargetDf.loc[:,'xpos'] = xlayout[i]
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            for j in range(0, tempDf.shape[0]):
                if tempDf.iloc[j,1] == 'other':
                    tempDf.iloc[j,1] = str(tempDf.iloc[j,0] + ' spp.')
            tempDf.columns = ['source','target','count']
            mean_val = [df[colour_col][df[cat_cols[i]]==x].mean() for x in df[cat_cols[i]]]
            tempDf.loc[:,'colour'] = mean_val
            tempDf.loc[:,'xpos'] = xlayout[i]
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
#         sourceTargetDf = sourceTargetDf.groupby(['source','target','colour','xpos']).agg({'count':'sum'}).reset_index()
        sourceTargetDf = sourceTargetDf.groupby(['source','target','xpos']).agg({'count':'sum', 'colour':'mean'}).reset_index()
        
    # add index for source-target pair
#     print(sourceTargetDf)
    labelList = list(set(sourceTargetDf.iloc[:,0:2].values.flatten()))
    labelList = list(dict.fromkeys(labelList))
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))
    
    return sourceTargetDf, colorList, labelList


# import taxonomy as tax
# taxa = tax.Taxonomy.from_ncbi(nodes_path='nodes.dmp', names_path='names.dmp')

def taxdist_old(query, suffix, reg_ids):
    file = 'blast/' + query + "." + suffix
    if checkfile(file) == 0:
        return
    blast = readblast('blast/' + query + "." + suffix)
    
    blast['species'] = ""
    blast['genus'] = ""
    blast['family'] = ""
    blast['order'] = ""
    blast['class'] = ""
    blast['phylum'] = ""
    blast['kingdom'] = ""
    blast['superkingdom'] = ""
    
    ranks = []
    
    for x in range(0, blast.shape[0]):
    #     print("Starting on ", x)
        taxid = blast['subject tax ids'][x]
        taxend = 0
        if str(taxid).find(";") != -1:
            taxid = str(taxid).split(";")[0]
        rank = ''
        if taxa.node(str(taxid)) == None:
            next
        else: 
            if taxa.node(str(taxid)).rank == 'no rank':
                taxid = taxa.node(str(taxid)).parent
            while taxend == 0:
                if taxa.node(str(taxid)) == None:
                    next
                else:
    #                 print("Processing taxids: ", taxid, rank, name)
                    name = taxa.node(str(taxid)).name
                    rank = taxa.node(str(taxid)).rank
                    if taxid in reg_ids:
                        blast.loc[x,'regulated'] = True
    #             print(rank)
                    if blast.columns.str.contains(rank).any():
                        blast.loc[x,rank] = name
                    taxid =  taxa.node(str(taxid)).parent
                    if name == "synthetic construct":
                        blast.loc[x,'species'] = " synthetic"
                        blast.loc[x,'genus'] = "synthetic"
                        blast.loc[x,'family'] = "synthetic "
                        blast.loc[x,'order'] = "synthetic  "
                        blast.loc[x,'class'] = "synthetic   "
                        blast.loc[x,'phylum'] = "synthetic    "
                        blast.loc[x,'kingdom'] = "synthetic     "
                        blast.loc[x,'superkingdom'] = "synthetic     "
                        taxend = 1
                        next
                    if rank == "superkingdom":
                        taxend = 1
                        next
    
    singletons = blast.species.value_counts().index[blast.species.value_counts()==1]
    blast['species_simplified'] = blast['species']
    blast.loc[blast.species.isin(singletons), 'species_simplified'] = 'other'
    
    return blast

def taxdist(query, suffix, reg_ids):
    file = 'blast/' + query + "." + suffix
    if checkfile(file) == 0:
        return
    blast = readblast('blast/' + query + "." + suffix)
    
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
    
    for x in range(0, blast.shape[0]):
        try:
            t = taxoniq.Taxon(blast['subject tax ids'][x])
            for level in t.ranked_lineage:
                if int(level.tax_id) in reg_ids:
                    blast.loc[x,'regulated'] = True
                if blast.columns.str.contains(level.rank.name).any():
                    blast.loc[x,level.rank.name] = level.scientific_name
                else:
                    blast[level.rank.name] = ""
                    blast.loc[x,level.rank.name] = level.scientific_name
        except:
            print('Taxon id', blast['subject tax ids'][x], 'is missing')
     
    singletons = blast.species.value_counts().index[blast.species.value_counts()==1]
    blast['species_simplified'] = blast['species']
    blast.loc[blast.species.isin(singletons), 'species_simplified'] = 'other'
    
    blast = blast.reset_index(drop=True)
    
    return blast

# REGULATED PATHOGEN LISTS

import gspread
from oauth2client.service_account import ServiceAccountCredentials
import taxoniq

def get_reg_ids():
    # use creds to create a client to interact with the Google Drive API
    scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name('../client-secret.json', scope)
    client = gspread.authorize(creds)
    
    reg_ids = []
    
    # read in bacteria
    sheet = client.open("Pathogen lists").worksheet("Bacteria")
    bact_list = sheet.col_values(1)
    australia = sheet.col_values(4)
    us_select = sheet.col_values(5)
    eu_export = sheet.col_values(7)
    russia = sheet.col_values(8)
    # check all the names are recognized
    # for bact in bact_list:
    #     t = taxoniq.Taxon(scientific_name=bact)
    #     print(t.scientific_name)
    
    # add taxon ids to regulated list
    for i in range(1,len(bact_list)):
        if ('Y' in [australia[i], us_select[i], eu_export[i], russia[i]]):
            t = taxoniq.Taxon(scientific_name=bact_list[i])
            reg_ids.append(t.tax_id)
    
    # read in viruses
    sheet = client.open("Pathogen lists").worksheet("Viruses")
    virus_list = sheet.col_values(1)
    australia = sheet.col_values(4)
    us_select = sheet.col_values(5)
    eu_export = sheet.col_values(7)
    russia = sheet.col_values(8)
    # check all the names are recognized
    # for vir in virus_list[1:]:
    #     t = taxoniq.Taxon(scientific_name=vir)
    #     print(t.scientific_name)
    
    # add taxon ids to regulated list
    for i in range(1,len(virus_list)):
        if ('Y' in [australia[i], us_select[i], eu_export[i], russia[i]]):
            t = taxoniq.Taxon(scientific_name=virus_list[i])
            reg_ids.append(t.tax_id)
            
    return reg_ids

def colourscale(reg_status, counts, averages):
    rmap = cm.get_cmap('OrRd', 100)
    nrmap = cm.get_cmap('Blues', 100)
    colours = []
    for i in range(len(reg_status)):
        colours.append('rgb' + str((np.array(nrmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')')) if reg_status[i] != counts[i] else colours.append('rgb' + str((np.array(rmap(averages[i]/averages.max()))*255)[:3].tolist()).replace('[', '(').replace(']', ')'))
    return colours
