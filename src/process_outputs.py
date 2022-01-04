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

reg_ids = pd.read_csv(os.environ['PFAMDB'] + '/biorisk/reg_taxids', header=None)


# taxonomic distribution plots
def plot_tax(query, reg_ids):
    blast = taxdist(query, reg_ids)
    if blast is not None:
        plot_pie(blast, query)
        # plot_sankey(blast, query, suffix)
        # plot_sunburst(blast, query, suffix)

# donut chart of taxon hits
def plot_pie(blast, query):
    levels = ['superkingdom', 'phylum', 'order', 'genus', 'species']
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
    fig.write_image("figures/" + query + "_tax_pie.png", width=700, height=700, scale=2)

##### plotting sequence alignments

# basic plotting of alignments of hits
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
        fig.add_shape(type="rect", name = hit_description, x0 = start, x1 = end, y0=0.5+i+1, y1=1.3+i+1, line=dict(color="white"), fillcolor=colour)
        fig.add_trace(go.Scatter(text = str(colour), x = [start, end, end, start, start], y=[0.5+i+1, 0.5+i+1, 1.3+i+1, 1.3+i+1, 0.5+i+1], fill="toself", line=dict(color="white")))
        fig.add_annotation(xanchor='left', text=hit_description, x=qlen, y=(0.5+i+1+1.3+i+1)/2,font=dict(family="Arial Narrow", size=int(100/yax), color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
    fig.update_xaxes(range=[0, qlen*2+3000])
    fig.update_yaxes(range=[yax+0.5, 0.5], showticklabels=False)
    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0))
    return fig

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
        fig.savefig(file + ".png")
        return
    
    hmmer = readhmmer(file)
    hmmer = hmmer.iloc[0:nhits,:]
    
    if hmmer.shape[0] < nhits:
        nhits = hmmer.shape[0]
    
    if re.search("biorisk", file):
        colours = colourscale([1.0] * hmmer.shape[0], [1.0] * hmmer.shape[0], pd.to_numeric(hmmer['score']))
    else:
        colours = colourscale([0.0] * hmmer.shape[0], [1.0] * hmmer.shape[0], pd.to_numeric(hmmer['score']))

    fig = plothits(hmmer["ali from"], hmmer["ali to"], hmmer['qlen'][0], hmmer["target name"], colours, nhits)
    fig.update_layout(showlegend=False, title={'text': 'HMMER Database Hits', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    fig.write_image("figures/" + file + ".png", width=1000, height=60*nhits+60, scale=2)

# plot BLAST results from --domtblout
def plot_blast(file, nhits=10):
    if checkfile(file) == 0:
        return
    if checkfile(file) == 2:
		# generate empty plot saying "no hits"
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig("figures/" + file + ".png")
        return
	
    if re.search(".nr.blastx", file):
        blast = taxdist(file, reg_ids)
    else:
        blast = readblast(file)
    blast = trimblast(blast)
    
    blast = blast.drop_duplicates('subject title') # drop hits with the same gene name
    blast = blast.reset_index()
    blast = blast.iloc[0:nhits,:]
    
    if blast.shape[0] < nhits:
        nhits = blast.shape[0]
    
    if re.search(".nr.blastx", file):
        colours = colourscale(blast['regulated'], [1.0] * blast.shape[0], pd.to_numeric(blast['% identity']))
#        print(colours)
    else:
        colours = colourscale([0.0] * blast.shape[0], [1.0] * blast.shape[0], pd.to_numeric(blast['% identity']))
#        print(colours)
    fig = plothits(blast['q. start'], blast['q. end'], blast['query length'][0], blast['subject title'], colours, nhits)
    fig.update_layout(showlegend=False, title={'text': 'BLAST Database Hits', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    fig.write_image("figures/" + file + ".png", width=1000, height=60*nhits+60, scale=2)




