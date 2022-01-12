import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as path_effects
from matplotlib import gridspec
import argparse
import pandas as pd
import numpy as np
import json
from utils import plot_style
import glob
import string
from collections import Counter

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--downsampleData')
	parser.add_argument('--config')
	parser.add_argument('--nIntroductions')
	args = parser.parse_args()
	plot_style()
	#args.downsampleData = 'data/19B_subclade/19B_subclade_sampled_5_rep_*_mcc_parsed.tsv'
	#args.config = 'config/downsampled_beast.json'
	plot_config = json.load(open(args.config, 'r'))
	# plot_config = json.load(open('config/downsampled_beast.json', 'r'))
	plot_config['effects'] = eval(plot_config['effects'])

	downsample_data = \
		pd.read_csv(args.downsampleData, sep='\t', header=None)
	plot_data = downsample_data.groupby(0).apply(lambda k: 
		k.sort_values(by=2).iloc[-5:,:]).reset_index(drop=True)
	plot_order = plot_data.groupby(1).apply(lambda k: 
			k[2].sum()).sort_values(
			ascending=False).index.values

	n_introductions = pd.read_csv(args.nIntroductions, sep='\t', header=None)	
	n_introductions = pd.read_csv('data/19B_subclade/19B_subclade_sampled_5_rep_*_location_tree_with_trait_importations_count.tsv', sep='\t', header=None)	
	
	n_introductions_counted = \
		{group_id: Counter(group[1]) for group_id, group in n_introductions.groupby(0)}
	n_introductions_counted = {group: 
		{key: val/sum(group_dat.values()) for 
			key, val in group_dat.items()} for group, group_dat in n_introductions_counted.items()}

	#fig, axs = plt.subplots(1, 3,figsize=(6.4*1.25, 4.8), constrained_layout=True, 
	#	gridspec_kw={'width_ratios': [0.1, 0.7, 0.2]})
	fig = plt.figure(figsize=(6.4*1.3, 4.8), constrained_layout=True)
	spec = fig.add_gridspec(ncols=11, nrows=5)

	ax0 = fig.add_subplot(spec[:,0])
	ax1 = fig.add_subplot(spec[:,1:8])
	

	axs = [ax0, ax1]
	for label_idx, label in enumerate(plot_order):
		axs[0].text(0, 0.2 + 0.075 * label_idx, plot_config['labels'][label],
				color=plot_config['colors'][label],
		        size=18, transform=axs[0].transAxes, 
		        path_effects=plot_config['effects'])


	axs[0].text(0, 0.2 + 0.075 * plot_order.shape[0], 'Other',
				color='#4d4d4d',
		        size=18, transform=axs[0].transAxes, 
		        path_effects=plot_config['effects'])


	_ = [axs[0].spines[loc].set_visible(False) for 
		loc in ['top', 'right', 'left', 'bottom']]

	axs[0].set_yticks([])
	axs[0].set_xticks([])

	ax=axs[1]
	print(plot_data[0].unique())
	ax.barh(plot_data[0].unique(), [1] * plot_data[0].unique().shape[0], 
		color='#4d4d4d', edgecolor='#333333', height=1)
	bottoms = [0] * plot_data[0].unique().shape[0]
	for loc in plot_order:
		loc_dat = plot_data[plot_data[1] == loc]
		loc_dat = loc_dat.set_index(0).reindex(plot_data[0].unique()).fillna(0)
		ax.barh(loc_dat.index, 
			loc_dat[2], 
			left=bottoms,
			facecolor=plot_config['colors'][loc], 
			edgecolor='#333333', height=1)
		bottoms += loc_dat[2].values


	ax.set_ylim(-0.5, plot_data[0].unique().shape[0]-0.5)
	ax.set_ylabel('Replicate')
	ax.set_xlabel('Posterior')
	ax.set_xlim(0, 1)

	group_axs = []
	for group_idx, (group_key, group_dat) in enumerate(n_introductions_counted.items()):
		group_ax = fig.add_subplot(spec[group_idx,8:])
		group_ax.bar(group_dat.keys(), group_dat.values(),
			facecolor=plot_config['colors']['GeorgiaUSA'],
			edgecolor='#333333', width=1.0)
		if group_idx == len(n_introductions_counted.keys())-1:
			group_ax.set_xlabel('Introductions')
		else:
			group_ax.set_xticklabels([])
		if group_idx == 2:
			group_ax.set_ylabel('Proportion')
		_ = [group_ax.spines[loc].set_visible(False) for loc in ['top', 'right']]
		group_axs.append(group_ax)
		#group_ax.set_aspect(1./group_ax.get_data_ratio())


	x_lims = [i.get_xlim() for i in group_axs]
	min_x = min([i[0] for i in x_lims])
	max_x = max([i[1] for i in x_lims])
	_ = [i.set_xlim(min_x, max_x) for i in group_axs]


	y_lims = [i.get_ylim() for i in group_axs]
	min_y = min([i[0] for i in y_lims])
	max_y = max([i[1] for i in y_lims])
	_ = [i.set_ylim(min_y, max_y) for i in group_axs]


	x_poss = [-0.12, -0.2]
	y_poss = [1.0, 1.0]
	for ax_idx, ax in enumerate(axs[1:]):
		ax.text(x_poss[ax_idx], y_poss[ax_idx], string.ascii_uppercase[ax_idx], 
			transform=ax.transAxes, 
			size=20, weight='bold', va="top")

	#fig.tight_layout()
	fig.savefig(plot_config['out_name'] + '.pdf')




if __name__ == "__main__":
    run()




