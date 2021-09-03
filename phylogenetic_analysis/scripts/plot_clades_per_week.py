import argparse
import baltic as bt
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import glob
import ast
from treetime.utils import numeric_date, datetime_from_numeric
import pandas as pd
import datetime
from collections.abc import Iterable
import json
import math
from utils import plot_style
import string
import numpy as np
from collections import Counter

def run():
parser = argparse.ArgumentParser()
parser.add_argument('--metadata')
parser.add_argument('--metadataDelim', default='\t')
parser.add_argument('--metadataCladeCol', default=2)
parser.add_argument('--metadataDateCol', default=1)
parser.add_argument('--config',
    help='json file with config information (labels, colors, etc.)')
args = parser.parse_args()
plot_style()
args.metadata = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted_GA_pangolin_dates.csv'
args.config='config/pango_lineage.json'

plot_config = json.load(open(args.config, 'r'))
plot_config['effects'] = eval(plot_config['effects'])


    metadata = pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
    print(metadata)
    metadata[args.metadataDateCol] = pd.to_datetime(metadata[args.metadataDateCol])
    metadata['week_end_date'] = pd.to_datetime(metadata[args.metadataDateCol]).apply(lambda k: k+datetime.timedelta(days= 6 - k.weekday()))
    clade_counts = Counter(metadata[args.metadataCladeCol])
    to_collapse = {key: f'Other ({key.split(".")[0]})' if value < plot_config['min_size'] else key for key, value in clade_counts.items() }
    metadata[args.metadataCladeCol] = metadata[args.metadataCladeCol].map(to_collapse)
n_per_week = metadata.groupby(['week_end_date', args.metadataCladeCol]).size().reset_index()
n_per_week = n_per_week.pivot(index='week_end_date', columns=2, values=0).fillna(0)
n_per_week.index = n_per_week.index.to_series().apply(numeric_date)
fig, ax = plt.subplots(constrained_layout=True, figsize=(6.4*1.25, 4.8))
cm = LinearSegmentedColormap.from_list('clades', [plot_config['colors'][i] for i in n_per_week.columns], N=n_per_week.shape[1])
    #(n_per_week.div(n_per_week.sum(axis=1), axis=0)*100).plot.area(ax=ax, colormap=cm, alpha=1.0, zorder=2)
n_per_week.plot.area(ax=ax, colormap=cm, alpha=1.0, zorder=2)
n_per_week['dt'] = n_per_week.index.to_series().apply(lambda k: datetime_from_numeric(k))
n_per_week.to_csv('.'.join(args.metadata.split('.')[0:-1])+'_n__per_week.csv')
x_ticks = [
    datetime.datetime.strptime('2020-03-01', '%Y-%m-%d'),
    datetime.datetime.strptime('2020-03-08', '%Y-%m-%d'),
    datetime.datetime.strptime('2020-03-15', '%Y-%m-%d'),
    datetime.datetime.strptime('2020-03-22', '%Y-%m-%d'),
    datetime.datetime.strptime('2020-03-31', '%Y-%m-%d')]
x_labels = \
    [i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
x_ticks = [numeric_date(i) for i in x_ticks]
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels)
ax.set_ylabel('Sequences/Wk', size=12)
ax.get_legend().remove()
ax.set_xlabel("Date (2020)")
x_pos = 0.8
y_pos = 0.9
ax_height = (ax.get_position().y1 - ax.get_position().y0) * ax.figure.bbox.height
        
for label_idx, label in enumerate(n_per_week.columns[:-1][::-1]):
    font_size = 14+int(ax.figure.bbox.height*0.01)*2.5
    text = ax.text(x_pos, y_pos - 0.075*label_idx, 
        label, color=plot_config['colors'][label],
                size=font_size, transform=ax.transAxes, 
                path_effects=plot_config['effects'])


fig.savefig(f'{plot_config["out_name"]}.pdf')




