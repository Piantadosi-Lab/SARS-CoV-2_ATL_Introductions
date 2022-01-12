#!/usr/bin/env python3
import argparse
import baltic as bt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import seaborn as sns
import glob
import ast
from treetime.utils import numeric_date, datetime_from_numeric
import pandas as pd
from datetime import datetime
from collections.abc import Iterable
import json
import math


def plot_style(grey='#333333'):
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['font.weight'] = 'light'
    mpl.rcParams['text.color'] = grey
    mpl.rcParams['axes.labelcolor'] = grey
    mpl.rcParams['xtick.color'] = grey
    mpl.rcParams['ytick.color'] = grey
    # Font sizes
    mpl.rcParams['figure.titlesize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    # Border colors
    mpl.rcParams['axes.edgecolor'] = grey


def import_tree(tree_file, plot_config):
	myTree = bt.loadNexus(tree_file, absoluteTime=False)
	# set any negative branch lengths to 0
	for k in myTree.Objects:
		if k.length < 0:
			k.length = 0
	myTree.drawTree()
	myTree.setAbsoluteTime(plot_config['max_date'])
	return(myTree)


def plot_divergence_tree(ax, tree_file=None, tree_name_sep='|', tree_name_field=1, trait_dict=None, plot_config=None, highlight_mrca=None):    
    myTree = bt.loadNewick(tree_file)
    myTree.drawTree()
    xMin = min([i.x for i in myTree.Objects])
    xMax = max([i.x for i in myTree.Objects])
    xSpan = xMax - xMin
    myTree.plotTree(ax, 
        x_attr=lambda k: k.x,
        width=plot_config['branch_width'] if 'branch_width' in plot_config.keys() else 1.0, 
        colour=plot_config['colors']['base'])
    ax.set_xlim((xMin-xSpan*0.1, xMax+xSpan*0.1))
    myTree.plotPoints(ax, 
        x_attr=lambda k: k.x, 
        size=lambda k: 40 if trait_dict[k.name.split(tree_name_sep)[tree_name_field]] in plot_config['colors'].keys() else 0.0, 
        colour=lambda k: plot_config['colors'][trait_dict[k.name.split(tree_name_sep)[tree_name_field]]] if 
            trait_dict[k.name.split(tree_name_sep)[tree_name_field]] in plot_config['colors'].keys() else 
            plot_config['colors']['base'], 
        outline_colour=plot_config['colors']['base'])
    # if we want to higlight the mrca, then get that MRCA
    if not (highlight_mrca is None):
    	mrca = myTree.commonAncestor(myTree.getExternal(lambda k: k.name.split(tree_name_sep)[tree_name_field] in highlight_mrca))
    	ax.scatter([mrca.x], [mrca.y], marker="P", zorder=5, 
    		edgecolor=plot_config['colors']['base'],
    		color=plot_config['colors']['highlight_mrca'],
    		s=500)
    ax.set_ylim(-myTree.ySpan*0.01, myTree.ySpan*1.01)
    ax.set_yticks([])
    ax.set_xticks([])
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left', 'bottom']]
    # add divergence scale
    # todo smarter and not manual
    ax.plot([xMin-xSpan*0.09, xMin-xSpan*0.09+0.0002], 
    	[myTree.ySpan*0.095, myTree.ySpan*0.095],
    	color=plot_config['colors']['base'])
    ax.plot([xMin-xSpan*0.09, xMin-xSpan*0.09], 
    	[myTree.ySpan*0.09, myTree.ySpan*0.1],
    	color=plot_config['colors']['base'])
    ax.plot([xMin-xSpan*0.09+0.0002, xMin-xSpan*0.09+0.0002], 
    	[myTree.ySpan*0.09, myTree.ySpan*0.1],
    	color=plot_config['colors']['base'])
    ax.text(xMin-xSpan*0.08, myTree.ySpan*0.08, '$2\\times10^{-4}$ subs/site', 
    	color=plot_config['colors']['base'], size=14, va='top')
    if 'labels' in plot_config.keys():
        ax_height = (ax.get_position().y1 - ax.get_position().y0) * ax.figure.bbox.height
        font_size = 10+int(ax.figure.bbox.height*0.01)*2.5
        starting_y = 0.13 + (1.2*font_size/ax_height)*(len(plot_config['labels'])+1)
        for label_idx, (label_key, label) in enumerate(plot_config['labels'].items()):
            text = ax.text(0, starting_y - (1.2*font_size/ax_height)*label_idx, 
                label, color=plot_config['colors'][label_key],
                        size=font_size, transform=ax.transAxes, 
                        path_effects=plot_config['effects'])
    return(ax)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--tree', 
	                    help='newick tree file to plot',
	                    default=None)
	parser.add_argument('--treeNameSep', 
	                    help='delimiter to seperate tree names on to match with metadata',
	                    default='|')
	parser.add_argument('--treeNameField', 
	                    help='which field in split tree names to match to metadata',
	                    default=1,
	                    type=int)
	parser.add_argument('--highlightMRCA')
	parser.add_argument('--metadata',
	    help='file with metadata')
	parser.add_argument('--metadataDelim',
	    help='file with metadata',
	    default='\t')
	parser.add_argument('--metadataRegionCol',
	    help='file with metadata',
	    type=int, 
	    default=7)
	parser.add_argument('--metadataIDCol',
	    help='column with id',
	    type=int, 
	    default=1)
	parser.add_argument('--config',
	    help='json file with config information (labels, colors, etc.)')
	args = parser.parse_args()
	plot_style()
	# reads in and processes config file
	plot_config = json.load(open(args.config, 'r'))
	plot_config['effects'] = eval(plot_config['effects'])
	plot_config['max_date'] = float(plot_config['max_date'])


	metadata = pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
	from collections import defaultdict
	trait_dict = defaultdict(lambda: None)
	trait_dict.update({i[args.metadataIDCol]: i[args.metadataRegionCol] for idx, i in metadata.iterrows()})

	print(trait_dict)
	highlight_mrca=None
	if args.highlightMRCA:
		highlight_mrca = \
			pd.read_csv(args.highlightMRCA, sep=args.metadataDelim,
				header=None)[args.metadataIDCol].values

	fig, ax = plt.subplots(1,1, figsize=(6.4, 4.8*2), constrained_layout=True)
	ax = plot_divergence_tree(ax, tree_file=args.tree, 
		tree_name_sep=args.treeNameSep, tree_name_field=args.treeNameField, 
		trait_dict=trait_dict, plot_config=plot_config, highlight_mrca=highlight_mrca)
	fig.savefig(f'{plot_config["out_name"]}.pdf')




if __name__ == "__main__":
    run()

