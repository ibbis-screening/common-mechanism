import pandas as pd
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import re
import os
from utils import *
import textwrap
import math

# reg_ids = pd.read_csv(args.database + '/biorisk/reg_taxids', header=None)
# vax_ids = pd.read_csv(args.database + '/benign/vax_taxids', header=None)

# taxonomic distribution plots
def plot_tax(file, reg_ids, query):
    blast = readblast(file)
    hits = taxdist(blast, reg_ids, query)
    if hits is not None:
        plot_pie(hits, query)

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
    fig.write_image(query + "_tax_pie.pdf", width=700, height=700, scale=2)

##### plotting sequence alignments

# basic plotting of alignments of hits
# def plothits_new(starts, ends, qlen, names, colours, nhits, max):
#     if nhits <5:
#         yax = 6
#     else:
#         yax=nhits+1
#     # fontsize = int(10/math.log(yax, 200))
#     xlim = qlen*2
#     fontsize = 20*(qlen/3000)
#     plts = make_subplots(rows=1, cols=2)
#     plts.add_trace(
#         fig = go.Figure(go.Scatter(x=[0,0], y=[0,0],mode="markers",marker=dict(color=[0,max], colorscale="Greys",colorbar=dict(title="Similarity", x=-0.15, xanchor="left"))), go.Layout(plot_bgcolor="white")),
#     row=1, col=1)
#     plts.add_trace(
#         fig.add_annotation(xanchor='left', text='Query', x=qlen, y=(0.5+1.3)/2,font=dict(family="Arial Narrow", size=fontsize, color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
#         for i in range(0, starts.shape[0]): # for each hit
#             fig.add_annotation(xanchor='left', text=hit_description, x=qlen, y=(0.5+i+1+1.3+i+1)/2,font=dict(family="Arial Narrow", size=fontsize, color="#000000"), bgcolor="#ffffff", showarrow=False), # Courier New, monospace,
#     row=1, col=2)
#     fig.add_shape(type="rect", name = 'Query', x0 = 1, x1 = qlen, y0=0.5, y1=1.3, line=dict(color="white"), fillcolor='grey')
#     for i in range(0, starts.shape[0]): # for each hit
#         start = starts[i]
#         end = ends[i]
#         max_name_length = 100
#         if len(names[i]) > max_name_length:
#             hit_description = names[i][:max_name_length] + "..."
#         else:
#             hit_description = names[i][:max_name_length]
#         hit_description = textwrap.fill(hit_description, 60).replace("\n", "<br>")
#         colour = colours[i]
#         fig.add_shape(type="rect", name = hit_description, x0 = start, x1 = end, y0=0.5+i+1, y1=1.3+i+1, line=dict(color="white"), fillcolor=colour)
#         fig.add_trace(go.Scatter(text = str(colour), x = [start, end, end, start, start], y=[0.5+i+1, 0.5+i+1, 1.3+i+1, 1.3+i+1, 0.5+i+1], fill="toself", line=dict(color="white")))
#         fig.add_annotation(xanchor='left', text=hit_description, x=qlen, y=(0.5+i+1+1.3+i+1)/2,font=dict(family="Arial Narrow", size=fontsize, color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
#     fig.update_xaxes(range=[0, xlim])
#     fig.update_yaxes(range=[yax+0.5, 0.5], showticklabels=False)
#     fig.update_layout(margin=dict(l=10, r=20, t=30, b=0))
#     return plts


# basic plotting of alignments of hits
def plothits(starts, ends, qlen, names, colours, nhits, max):
    if nhits <5:
        yax = 6
    else:
        yax=nhits+1
    fig = go.Figure(go.Scatter(x=[0,0], y=[0,0],mode="markers",marker=dict(color=[0,max], colorscale="Greys",colorbar=dict(title="Similarity", x=-0.15, xanchor="left"))), go.Layout(plot_bgcolor="white"))
    fig.add_shape(type="rect", name = 'Query', x0 = 1, x1 = qlen, y0=0.5, y1=1.3, line=dict(color="white"), fillcolor='grey')
    fig.add_annotation(xanchor='left', text='Query', x=qlen, y=(0.5+1.3)/2,font=dict(family="Arial Narrow", size=int(100/yax), color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
    for i in range(0, starts.shape[0]): # for each hit
        start = starts[i]
        end = ends[i]
        max_name_length = 100
        if len(names[i]) > max_name_length:
            hit_description = names[i][:max_name_length] + "..."
        else:
            hit_description = names[i][:max_name_length]
        colour = colours[i]
        fig.add_shape(type="rect", name = hit_description, x0 = start, x1 = end, y0=0.5+i+1, y1=1.3+i+1, line=dict(color="white"), fillcolor=colour)
        fig.add_trace(go.Scatter(text = str(colour), x = [start, end, end, start, start], y=[0.5+i+1, 0.5+i+1, 1.3+i+1, 1.3+i+1, 0.5+i+1], fill="toself", line=dict(color="white")))
        fig.add_annotation(xanchor='left', text=hit_description, x=qlen, y=(0.5+i+1+1.3+i+1)/2,font=dict(family="Arial Narrow", size=int(100/yax), color="#000000"), bgcolor="#ffffff", showarrow=False) # Courier New, monospace
    fig.update_xaxes(range=[0, qlen*2+1000])
    fig.update_yaxes(range=[yax+0.5, 0.5], showticklabels=False)
    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0))
    return fig


# plot HMMER results from --domtblout
def plot_hmmer(file, lookup, nhits=10):
    if check_blastfile(file) == 0:
        return
    if check_blastfile(file) == 2:
		# generate empty plot saying "no hits"
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5, 0.5, 'No hits', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig(os.path.abspath(file + ".png"))
        return
    
    hmmer = readhmmer(file)
    hmmer = trimhmmer(hmmer)
    hmmer = hmmer.iloc[0:nhits,:]
     
    # check if this is a batch of sequences
    seq_names = hmmer['query name']
    seq_pres = set([i.split('_')[0] for i in seq_names])
  
    if hmmer.shape[0] < nhits:
        nhits = hmmer.shape[0]
    
    if re.search("biorisk", file):
        colours = colourscale([1.0] * hmmer.shape[0], [1.0] * hmmer.shape[0], pd.to_numeric(hmmer['score']))
    else:
        colours = colourscale([0.0] * hmmer.shape[0], [1.0] * hmmer.shape[0], pd.to_numeric(hmmer['score']))

    new_names = []
    for model in range(hmmer.shape[0]):
        name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
        # print(name_index)
        # hmmer['description'][model] = lookup['Description'][name_index[0]]
        try:
            new_names.append(lookup['Description'][name_index[0]])
        except:
            new_names.append("")
        # print(lookup['Description'][name_index[0]])
    hmmer['description'] = new_names
    
    # print(hmmer)
    fig = plothits(hmmer["ali from"], hmmer["ali to"], hmmer['qlen'][0], hmmer["description"], colours, nhits, hmmer['score'].max())
    fig.update_layout(showlegend=False, title={'text': 'HMMER Database Hits', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    # print("Path specified: ", file + ".png")
    fig.write_image(file + ".png", width=1000, height=60*nhits+60, scale=2)

# plot BLAST results
def plot_blast(file, reg_ids, vax_ids, nhits):
    if check_blastfile(file) == 0:
        return
    if check_blastfile(file) == 2:
		# generate empty plot saying "no hits"
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig(os.path.abspath(file + ".png"))
        return
	
    blast = readblast(file)
    blast = taxdist(blast, reg_ids, vax_ids)
    
    blast = blast.sort_values(by=['% identity', 'bit score'], ascending=False)
    print(blast[['subject acc.', 'subject tax ids', 'regulated', '% identity', 'q. start', 'q. end', 's. start', 's. end']])
    blast = trimblast(blast)
    print(blast[['subject acc.', 'subject tax ids', '% identity', 'q. start', 'q. end', 's. start', 's. end']])
    blast = tophits(blast)
    print(blast[['subject acc.', 'subject tax ids', 'regulated', '% identity', 'q. start', 'q. end', 's. start', 's. end']])
    
    blast = blast.drop_duplicates('subject title') # drop hits with the same gene name
    blast = blast.reset_index()
    blast = blast.iloc[0:nhits,:]
        
    if blast.shape[0] < nhits:
        nhits = blast.shape[0]
    
    colours = colourscale(blast['regulated'], [1.0] * blast.shape[0], pd.to_numeric(blast['% identity']))
    names = blast['subject acc.'] + ": " + blast['subject title']
    fig = plothits(blast['q. start'], blast['q. end'], blast['query length'][0], names, colours, nhits, 100)
    fig.update_layout(showlegend=False) # , title={'text': 'Database Hits', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'}
    fig.write_image(os.path.abspath(file + ".png"), width=1000, height=60*nhits+90, scale=2)





