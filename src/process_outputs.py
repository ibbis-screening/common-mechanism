import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import plotly
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib
import re
import os
import subprocess
from Bio import ExPASy
from Bio import SwissProt
from utils import *
# pio.renderers.default = "browser"

lookup = pd.read_csv("../databases/VFDB_setA_geneids", sep="\t", header=None)
lookup.head()

# taxonomic distribution plots
def plot_tax(query, reg_ids):
    blast = taxdist(query, suffix, reg_ids)
    if blast is not None:
        plot_pie(blast, query, suffix)
        # plot_sankey(blast, query, suffix)
        # plot_sunburst(blast, query, suffix)

# donut chart of taxon hits
def plot_pie(blast, query, suffix):
    levels = ['superkingdom', 'phylum', 'class', 'genus', 'species']
    traces = []
    sizes = [1,3.5,9,16,25,34]
    holes = [0,0.5,0.62,0.74,0.8,0.85]
    for i in range(len(levels)): 
        blast.loc[blast[levels[i]]=='',levels[i]] = blast.loc[blast[levels[i]]=='',levels[i-1]] + "  "
        blast_collapse = blast.sort_values(levels).groupby(levels[:(i+1)], sort=False)
        averages = blast_collapse.mean()['% identity']
        counts = list(blast_collapse.count()['% identity'])
        colours = colourscale(blast_collapse.sum()['regulated'], counts, averages)
        values = list(blast_collapse.count()['% identity'])
        values_scaled = np.asarray(values)/sum(np.asarray(values))*sizes[i]
        labels = np.concatenate(list(blast_collapse[levels[i]].unique()))
        trace = go.Pie(
            hole = holes[i],
            scalegroup = 'one',
            sort = False,
            direction='clockwise',
            values = values_scaled,
            labels = labels,
            marker={'colors': colours}
        )
        traces.append(trace)
    fig = go.Figure(data=traces)
    fig.update_traces(textposition='inside', textinfo='label')
    fig.update_layout(showlegend=False)
    fig.write_image("figures/" + query + "." + suffix + "_tax_pie.png", width=700, height=700, scale=2)

# Sankey diagram of taxon hits
def plot_sankey(blast, query, suffix):
    blast['count'] = 1
    sourceTargetDf, colorList, labelList = genSankey(blast, cat_cols=['superkingdom', 'phylum', 'class', 'genus', 'species_simplified'], value_cols='count', colour_col='log evalue', title='Taxonomic distribution of hits')
    
    norm = matplotlib.colors.Normalize(vmin=10, vmax=300)
    cmap = cm.YlOrRd
    m = cm.ScalarMappable(cmap=cmap)
    
    sourceTargetDf = sourceTargetDf[sourceTargetDf['count']>3]
    
    fig = go.Figure(data=[go.Sankey(
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(color = "black", width = 0.5),
          label = labelList,
          color = "white" # colorList # colours of the nodes
        ),
        link = dict(
            source = sourceTargetDf['sourceID'],
            target = sourceTargetDf['targetID'],
            value = sourceTargetDf['count'],
            color = [matplotlib.colors.rgb2hex(x) for x in m.to_rgba(sourceTargetDf['colour'])]
      ))])
    # plotly.offline.plot(fig, validate=False)
    fig.update_layout(title_text="Taxonomic distribution of hits", font_size=20, autosize=False, width=1500, height=700)
    fig.write_image("figures/" + query + "." + suffix + "_tax_sankey.png", width=1500, height=700, scale=2)

# plot HMMER results from --domtblout
def plot_hmmer(file, nhits=10):
    if checkfile(file) == 0:
        return
    if checkfile(file) == 2:
		# generate empty plot saying "no hits"
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig("file.png")
        return
	
	hmmer = read_hmmer(file)

# processing VFDB hits
def vfdb(query, nhits=20):
    # file = "blast/" + query + ".vfdb.blastp"
    file = "blast/" + query + ".vfdb.blastx"
    if checkfile(file) == 0:
        print('Zero status')
        return
    if checkfile(file) == 2:
        print("No hits")
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig("figures/" + query + "_vfdb.png")
        return
    
    blast = readblast("blast/" + query + ".vfdb.blastx")
    blast = blast.iloc[0:nhits,:]
    blast['name'] = [re.sub(r"VFG\d+\(gb\|\S+\)\ ","",x) for x in blast['subject title']]

    funcs = pd.read_csv("../databases/VFs.csv")
    funcs['include'] = False
    for i in range(funcs.shape[0]):
        try:
            if taxoniq.Taxon(scientific_name=funcs.loc[i,'Bacteria']).tax_id in reg_ids:
                funcs.loc[i,'include'] = True
        except:
            'Whatever'
    funcs = funcs[funcs['include'] == True]
    
    # visualising hit lengths and homologies
    hit_descriptions = []
    for i in range(0, blast.shape[0]):
        func_an = blast['subject title'][i]
        if (re.search("\(VF\d+\)", func_an)) != None:
            result = re.search("\(VF\d+\)", func_an)
            vf = result.group(0)
            functions = funcs['Keyword'][funcs['VFID'].str.contains(vf)].to_string(index=False)
            short_name = funcs['VF_Name'][funcs['VFID'].str.contains(vf)].to_string(index=False)
            hit_descriptions.append(str(short_name + " : " + functions))
        else:
            hit_descriptions.append(func_an)
    
    fig = plothits(blast["q. start"], blast["q. end"], blast['query length'][0], blast["subject title"], pd.to_numeric(blast['% identity']), nhits)
    fig.update_layout(showlegend=False, title={'text': 'Virulence Factor Database functional annotation', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    fig.write_image("figures/" + query + "_vfdb.png", width=1000, height=60*nhits+60, scale=2)
    
    return blast

# processing blacklist hits
def biorisk(query, nhits=5):
    file = "blast/" + query + ".biorisk.tblastx"
    if checkfile(file) == 0:
        return
    if checkfile(file) == 2:
        print("No hits")
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits to biorisk database', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig("figures/" + query + "_biorisk.png")
        return
    blast = readblast("blast/" + query + ".biorisk.tblastx")
    blast = trimblast(blast)
    blast = blast.iloc[0:nhits,:]
    blast['flag'] = sum([blast['% identity']>0.8, blast['evalue']<0.0001, blast['q. coverage']>0.8, blast['s. coverage']>0.8])
    # nhits = blast.shape[0]+1
    
    # visualising hit lengths and homologies
    fig = plothits(blast["q. start"], blast["q. end"], blast['query length'][0], blast["subject title"], pd.to_numeric(blast['% identity']), nhits)
    fig.update_layout(title={'text': 'Biorisk database hits', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'}, showlegend=False)
    fig.write_image("figures/" + query + "_biorisk.png", width=1000, height=60*nhits+60, scale=2) # blast.shape[0]*20+200
    
    return blast

# processing Swissprot results
def sprot(query, nhits=5):
    file = 'blast/' + query + ".sprot.blastx"
    if checkfile(file) == 0:
        return
    if checkfile(file) == 2:
        print("No hits")
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits to SwissProt', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig("figures/" + query + "_sprot.png")
        return
    blast = readblast('blast/' + query + ".sprot.blastx")
    blast = trimblast(blast)
    blast = blast.iloc[0:nhits,:]
    names = []
    for x in blast['subject title']:
        if (re.search(r"RecName: Full=.*(;|\[)", x)) != None:
            result = re.search(r"RecName: Full=.*(;|\[)", x)
            name = result.group(0)
            name = re.sub("\[", "", name)
            name = re.sub("RecName: Full=", "", name)
            name = re.sub("; Short.*", "", name)
            names.append(name)
        else:
            print("Name error: ", x)
    
    blast['flag'] = ''
    funsocs = []
    goterms = []
    rel_keys = []
    oth_keys = []
    
    # nhits = blast.shape[0]+1
    
    for i in range(0, blast.shape[0]):
        sprotid = blast.iloc[i,2]
        handle = ExPASy.get_sprot_raw(sprotid)
        record = SwissProt.read(handle)
        rel_keys = add_rel_keys(rel_keys, record)
        oth_keys = add_keys(oth_keys, record)
        if rel_keys[i] != '':
            blast.loc[i,'flag'] = True
        else:
            blast.loc[i,'flag'] = False

    # print('keywords: ', keys)
    desc_format = "{:<8} {:<20} {:<20}"
    descriptions = []
    for i in range(blast.shape[0]):
        if rel_keys[i] != '':
            print('Virulence function')
        sep = ' | '
        # print(i)
        # descriptions.append(desc_format.format(blast.iloc[i,2], funsocs[i], goterms[i], names[i]))
        # descriptions.append(desc_format.format(blast.iloc[i,2], rel_keys[i], oth_keys[i]))
        descriptions.append(sep.join([blast.iloc[i,2], rel_keys[i], names[i]]))
    if len(descriptions) < blast.shape[0]:
        print('Error - length of descriptions does not match length of hits', descriptions)
    # descriptions = [i + " " + j for i, j in zip(funsocs, blast['subject title'])]
    
    fig = plothits(blast['q. start'], blast['q. end'], blast['query length'][0], descriptions, pd.to_numeric(blast['% identity']), nhits)
    fig.update_layout(showlegend=False, title={'text': 'Swissprot hits and relevant functional annotations', 'y':1, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    fig.write_image("figures/" + query + "_sprotfunsocs.png", width=1000, height=60*nhits+60, scale=2)
    
    return blast

# TrEMBL results
def trembl(query, nhits=5):
    file = 'blast/' + query + ".trembl.blastx.pool"
    if checkfile(file) == 0:
        return
    if checkfile(file) == 2:
        print("No hits")
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits to TrEMBL', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig("figures/" + query + "_trembl.png")
        return
    blast = readblast(file)
    blast = trimblast(blast)
    blast = blast.iloc[0:nhits,:]
    names = []
    
    blast['flag'] = ''
    funsocs = []
    goterms = []
    rel_keys = []
    oth_keys = []
    
    # nhits = blast.shape[0]+1
    
    for i in range(0, blast.shape[0]):
        sprotid = blast.iloc[i,2]
        handle = ExPASy.get_sprot_raw(sprotid)
        record = SwissProt.read(handle)
        rel_keys = add_rel_keys(rel_keys, record)
        oth_keys = add_keys(oth_keys, record)
        if rel_keys[i] != '':
            blast.loc[i,'flag'] = True
        else:
            blast.loc[i,'flag'] = False

    # print('keywords: ', keys)
    desc_format = "{:<8} {:<20} {:<20}"
    descriptions = []
    for i in range(blast.shape[0]):
        if rel_keys[i] != '':
            print('Virulence function')
        sep = ' | '
        # print(i)
        # descriptions.append(desc_format.format(blast.iloc[i,2], funsocs[i], goterms[i], names[i]))
        # descriptions.append(desc_format.format(blast.iloc[i,2], rel_keys[i], oth_keys[i]))
        descriptions.append(sep.join([blast.iloc[i,2], rel_keys[i], names[i]]))
    if len(descriptions) < blast.shape[0]:
        print('Error - length of descriptions does not match length of hits', descriptions)
    # descriptions = [i + " " + j for i, j in zip(funsocs, blast['subject title'])]
    
    fig = plothits(blast['q. start'], blast['q. end'], blast['query length'][0], descriptions, pd.to_numeric(blast['% identity']), nhits)
    fig.update_layout(showlegend=False, title={'text': 'Swissprot hits and relevant functional annotations', 'y':1, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    fig.write_image("figures/" + query + "_sprotfunsocs.png", width=1000, height=60*nhits+60, scale=2)
    
    return blast

# functions to add annotations to BLAST hit plots
def add_keys(keylist, record):
    keywords = record.keywords
    f=open("keywords.txt", "a+")
    sep = ', '
    f.write(sep.join(keywords) + "\n")
    r = re.compile('^(?!Toxin|Virulence|DNA replication)')
    relevant_keys = list(filter(r.match, keywords))
    keylist.append(sep.join(relevant_keys))
    return keylist

def add_rel_keys(keylist, record):
    keywords = record.keywords
    r = re.compile('Toxin|Virulence|DNA replication')
    relevant_keys = list(filter(r.match, keywords))
    sep = ', '
    keylist.append(sep.join(relevant_keys))
    return keylist

def add_go(golist, record):
    goterms = []
    gos = ""
    for ex_db_data in record.cross_references:
        extdb,extid,goterm = ex_db_data[:3]
        if extdb=="GO":
            if re.match("P:", goterm):
                gos = goterm + ", " + gos
                goterms.append(gos)
    sep = ','
    golist.append(sep.join(gos))
    return golist

def add_funsocs(funlist, sprotid):
    funsocs = []
    cmd = "grep " + sprotid + " /Users/nw17/Documents/NTI/databases/funsocs.tsv"
    if subprocess.call(cmd, shell=True) == 0:
        returned_output = subprocess.check_output(cmd, shell = True)
        hit = returned_output.decode("utf-8").split(sep="\t")
        if hit.count("1")>0:
            indices = [hit.index(i) for i in ["1"]]
            funs = ['Entry','nonviral_adhesion','viral_adhesion','antibiotic_resistance','avirulence_plant','counter_immunoglobulin','cytotoxicity','develop_in_host','disable_organ','plant_rna_silencing','induce_inflammation','nonviral_invasion','viral_invasion','degrade_ecm','host_cytoskeleton','host_cell_cycle','host_cell_death','host_gtpase','host_transcription','host_translation','host_ubiquitin','host_xenophagy','nonviral_movement','viral_movement','toxin_synthase','resist_complement','resist_oxidative','secreted_effector','secretion','suppress_detection','bacterial_counter_signaling','viral_counter_signaling','virulence_regulator','virulence_activity','point_total'][indices[0]]
            funsocs.append(funs)
        else: 
            funsocs.append("No funsocs")
    else:
        funsocs.append("No funsocs")
    sep = ','
    funlist.append(sep.join(funsocs))
    return funlist

# read in queries mapped to sequence database
def readmap(query):
    mapping = pd.read_csv("../queries/" + query + ".trim", sep=" ", comment='@', header=None)
    mapped = mapping[mapping[2]!="*"]
    for ref in set(mapped[2]):
        hit = mapped[mapped[2]==ref]
        name = [i for i in lookup[0] if ref in i][0]
        fig = go.Figure(go.Scatter(x=[0,0], y=[0,0], fill="none", mode="markers", marker=dict(color="white")), go.Layout(plot_bgcolor="white"))
        for i in range(0, hit.shape[0]):
            start = int(hit[3].iloc[i])
            matchinfo = re.findall(r"\d+\w", str(hit.iloc[i,5]))
            for sect in matchinfo:
                if re.findall(r"[A-Z]", sect)[0] == "S":
                    colour="grey"
                if re.findall(r"[A-Z]", sect)[0] == "M":
                    colour="blue"
                    length = re.findall(r"\d+", sect)[0]
                    end = int(start) + int(length)
                    fig.add_shape(type="rect", name = "", x0 = start, x1 = end, y0=0.5+i, y1=1.3+i, line=dict(color=None), fillcolor=colour)
                    fig.add_trace(go.Scatter(text = name, x = [start, end, end, start, start], y=[0.5+i, 0.5+i, 1.3+i, 1.3+i, 0.5+i], fill="toself", line=dict(color="white")))
                    start = end
        fig.update_layout(xaxis_title=name, showlegend=False, margin=dict(l=0, r=0, t=0, b=0),font=dict(size=10))
        fig.write_image("figures/" + query + "_map.png", width=500, height=200, scale=2)

# create a summary figure of pipeline results
from PIL import Image
import cv2
import os

def results_page(query):
    imgs = []
    cvs = Image.new('RGB',(3000,1000), (255, 255, 255))
    # taxonomic distribution
    if os.path.isfile('figures/' + query + '.nr.blastx_tax_pie.png'):
        im = Image.open('figures/' + query + '.nr.blastx_tax_pie.png')
        im = im.resize((900,900))
        cvs.paste(im, (0,0))
    else:
        print("No taxonomic distribution image")
    # blacklist
    if os.path.isfile('figures/' + query + '_blacklist.png'):
        im = Image.open('figures/' + query + '_blacklist.png')
        im = im.resize((2000,500))
        cvs.paste(im, (1000,0))
    else:
        print("No blacklist image")
    # functional annotation
    if os.path.isfile('figures/' + query + '_sprotfunsocs.png'):
        im = Image.open('figures/' + query + '_sprotfunsocs.png')
        im = im.resize((2000,500))
        cvs.paste(im, (1000,500))
    else:
        print("No annotation image")
    return cvs
