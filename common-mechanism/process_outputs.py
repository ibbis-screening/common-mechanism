import os
import re
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
from utils import colourscale, is_empty, has_hits,readblast, readhmmer, taxdist, tophits, trimblast, trimhmmer

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
    for i, level in enumerate(levels):
        blast.loc[blast[level]=='',level] = blast.loc[blast[level]=='',levels[i-1]] + "  "
        blast_collapse = blast.sort_values(levels).groupby(levels[:(i+1)], sort=False)
        averages = blast_collapse['% identity'].mean()
        counts = list(blast_collapse.count()['% identity'])
        colours = colourscale(blast_collapse['regulated'].sum(), counts, averages)
        values = list(blast_collapse.count()['% identity'])
        values_scaled = np.asarray(values)/sum(np.asarray(values))*sizes[i]
        labels = np.concatenate(list(blast_collapse[level].unique()))
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
def plothits(starts, ends, qlen, names, colours, nhits, max):
    if nhits <5:
        yax = 6
    else:
        yax=nhits+1
    fig = go.Figure(go.Scatter(x=[0,0],
                               y=[0,0],
                               mode="markers",
                               marker={"color": [0,max]},
                               colorscale="Greys",
                               colorbar={"title": "Similarity", "x": -0.15, "xanchor":"left"}),
                    go.Layout(plot_bgcolor="white"))
    fig.add_shape(type="rect",
                  name = 'Query',
                  x0 = 1,
                  x1 = qlen,
                  y0 = 0.5,
                  y1 = 1.3,
                  line = {"color": "white"},
                  fillcolor = 'grey')
    fig.add_annotation(xanchor = 'left',
                       text='Query',
                       x=qlen,
                       y=(0.5+1.3)/2,
                       font={"family": "Arial Narrow", "size": int(100/yax), "color": "#000000"},
                       bgcolor="#ffffff",
                       showarrow=False)

    for i in range(0, starts.shape[0]): # for each hit
        start = starts[i]
        end = ends[i]
        max_name_length = 600
        if len(names[i]) > max_name_length:
            hit_description = names[i][:max_name_length] + "..."
        else:
            hit_description = names[i][:max_name_length]
        # hit_description = names[i]
        colour = colours[i]
        fig.add_shape(type="rect",
                      name = hit_description,
                      x0 = start,
                      x1 = end,
                      y0=0.5+i+1,
                      y1=1.3+i+1,
                      line={"color": "white"},
                      fillcolor=colour)
        fig.add_trace(go.Scatter(text = str(colour),
                                 x = [start, end, end, start, start],
                                 y = [0.5+i+1, 0.5+i+1, 1.3+i+1, 1.3+i+1, 0.5+i+1],
                                 fill="toself",
                                 line= {"color" :"white"}))
        fig.add_annotation(xanchor='left',
                           text=hit_description,
                           x=qlen,
                           y=(0.5+i+1+1.3+i+1)/2,
                           font={"family": "Arial Narrow", "size": (100/yax), "color": "#000000"},
                           bgcolor="#ffffff",
                           showarrow=False)
    fig.update_xaxes(range=[0, qlen*2+1000])
    fig.update_yaxes(range=[yax+0.5, 0.5], showticklabels=False)
    fig.update_layout(margin=dict(l=10, r=20, t=30, b=0))
    return fig


# plot HMMER results from --domtblout
def plot_hmmer(file, lookup, nhits=10):
    if is_empty(file):
        return
    if not has_hits(file):
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
    if hmmer.shape[0] < nhits:
        nhits = hmmer.shape[0]

    if re.search("biorisk", file):
        colours = colourscale([1.0] * hmmer.shape[0], [1.0] * hmmer.shape[0], pd.to_numeric(hmmer['score']))
    else:
        colours = colourscale([0.0] * hmmer.shape[0], [1.0] * hmmer.shape[0], pd.to_numeric(hmmer['score']))

    new_names = []
    for model in hmmer['target name']:
        new_names.append(lookup['Description'][lookup['ID'] == model].iloc[0])
    hmmer['description'] = new_names

    fig = plothits(hmmer["ali from"],
                   hmmer["ali to"],
                   hmmer['qlen'][0],
                   hmmer["description"],
                   colours,
                   nhits,
                   hmmer['score'].max())
    fig.update_layout(showlegend=False,
                      title={'text': 'HMM Database Hits', 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'})
    fig.write_image(file + ".png", width=1000, height=60*nhits+60, scale=2)

# plot BLAST results
def plot_blast(file, reg_ids, vax_ids, nhits):
    if is_empty(file):
        return
    if not has_hits(file):
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

    blast = trimblast(blast)
    blast = tophits(blast)

    blast = blast.drop_duplicates('subject acc.') # drop hits with the same gene name
    blast = blast.reset_index()

    blast['label'] = blast.groupby(['q. start','q. end'])['subject acc.'].transform(lambda x: ','.join(x))
    blast['description'] = blast.groupby(['q. start','q. end'])['subject title'].transform(lambda x: ','.join(x))
    blast['regulated_sum'] = blast.groupby(['q. start','q. end'])['regulated'].transform(lambda x: any(x))
    blast = blast.drop_duplicates(['q. start','q. end'])
    blast = blast.sort_values(by=['q. start'], ascending=True)
    blast.reset_index(inplace=True)

    if blast.shape[0] < nhits:
        nhits = blast.shape[0]

    colours = colourscale(blast['regulated_sum'], [1.0] * blast.shape[0], pd.to_numeric(blast['% identity']))
    names = blast['description']

    blast = blast.iloc[0:nhits,:]
    names = names[0:(nhits)]
    colours = colours[0:(nhits)]

    fig = plothits(blast['q. start'], blast['q. end'], blast['query length'][0], names, colours, nhits, 100)
    fig.update_layout(showlegend=False)
    fig.write_image(os.path.abspath(file + ".png"), width=1000, height=60*nhits+90, scale=2)



# plot BLAST results from fragmented noncoding file
def plot_blast_frag(file, reg_ids, vax_ids, nhits):
    if is_empty(file):
        return
    if not has_hits(file):
		# generate empty plot saying "no hits"
        fig = plt.figure(figsize=(10,3))
        ax = fig.add_subplot(111, frameon=False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.text(0.5,0.5, 'No hits', fontsize=30, verticalalignment='center', horizontalalignment='center')
        fig.savefig(os.path.abspath(file + ".png"))
        return

    blast = readblast(file)
    blast['coords'] = blast['query acc.'].transform(lambda x: re.search(':(.+?)$', x).group(1))
    coords = blast["coords"].str.split("-", n = 1, expand = True)
    blast['start_add'] = coords[0]
    blast['start_add'] = blast['start_add'].transform(lambda x: int(x))
    blast['q. start'] = blast['q. start'] + blast['start_add']
    blast['q. end'] = blast['q. end'] + blast['start_add']
    blast['query length'] = blast['q. end'].max()

    blast = taxdist(blast, reg_ids, vax_ids)

    blast = blast.sort_values(by=['% identity', 'bit score'], ascending=False)
    blast = trimblast(blast)
    blast = tophits(blast)

    blast = blast.drop_duplicates('subject acc.') # drop hits with the same gene name
    blast = blast.reset_index()

    blast['label'] = blast.groupby(['q. start','q. end'])['subject acc.'].transform(lambda x: ','.join(x))
    blast['description'] = blast.groupby(['q. start','q. end'])['subject title'].transform(lambda x: ','.join(x))
    blast['regulated_sum'] = blast.groupby(['q. start','q. end'])['regulated'].transform(lambda x: any(x))
    blast = blast.drop_duplicates(['q. start','q. end'])
    blast.reset_index(inplace=True)

    if blast.shape[0] < nhits:
        nhits = blast.shape[0]

    colours = colourscale(blast['regulated_sum'], [1.0] * blast.shape[0], pd.to_numeric(blast['% identity']))
    names = blast['description']

    blast = blast.iloc[0:nhits,:]
    names = names[0:(nhits)]
    colours = colours[0:(nhits)]

    fig = plothits(blast['q. start'], blast['q. end'], blast['query length'][0], names, colours, nhits, 100)
    fig.update_layout(showlegend=False) 
    fig.write_image(os.path.abspath(file + ".png"), width=1000, height=60*nhits+90, scale=2)
