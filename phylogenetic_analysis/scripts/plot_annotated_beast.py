import argparse
import baltic as bt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from treetime.utils import numeric_date
from matplotlib.colors import LinearSegmentedColormap
import string
import json
from datetime import datetime
import pandas as pd
from utils import plot_style


def import_tree(tree_file, plot_config):
	myTree = bt.loadNexus(tree_file, absoluteTime=False)
	prev_x_locs = {i.name: i.x for i in myTree.Objects if i.branchType=='leaf'}
	# set any negative branch lengths to 0
	for k in myTree.Objects:
		if k.length < 0:
			#delta = -1*k.length
			for child in k.children:
				child.length = child.length + k.length
			k.length = 0
	new_x_locs = {i.name: i.x for i in myTree.Objects if i.branchType=='leaf'}
	for key in prev_x_locs.keys():
		if round(prev_x_locs[key], 8) != round(new_x_locs[key], 8):
			raise Exception('tip locations changed')
	myTree.drawTree()
	# collapse negative branches
	myTree.setAbsoluteTime(plot_config['max_date'])
	#myTree.collapseBranches(collapseIf=lambda x:x.length<=0.0,verbose=True)
	return(myTree)



def get_lims(myTree):
	xMin = min([i.absoluteTime for i in myTree.Objects])
	xMax = max([i.absoluteTime for i in myTree.Objects])
	xSpan = xMax - xMin
	yMin = min([i.y for i in myTree.Objects])
	yMax = max([i.y for i in myTree.Objects])
	ySpan = yMax - yMin
	return(xSpan, (xMin, xMax), ySpan, (yMin, yMax))


def annotate_nodes(ax, mySubTree, xSpan, ySpan, nodes_to_annotate, plot_config):
	for node in nodes_to_annotate: 
		print(node)
		# plot node timing confidence intervals
		print(plot_config['max_date'] - node.traits["height"])
		xmin = plot_config['max_date'] - node.traits["height_95%_HPD"][1]
		xmax = plot_config['max_date'] - node.traits["height_95%_HPD"][0]
		print(xmin)
		print(xmax)
		ax.add_patch(Rectangle((xmin,node.y-0.01*mySubTree.ySpan),
			xmax-xmin, 0.02*mySubTree.ySpan,
			facecolor=plot_config['colors'][node.traits['location']],
			edgecolor=plot_config['colors'][node.traits['location']], 
			alpha=0.5)) 
		# plot pie chart showing reconstruction of node states
		location_prob = \
			[i for i in zip(node.traits['location.set'], 
				node.traits['location.set.prob']) if 
				i[0] in plot_config['labels'].keys()]
		# can assign all "other" to USA because we have colors for all other international locations
		location_prob.append(('Other (USA)', 1-sum([i[1] for i in location_prob])))
		print(location_prob)
		print(sum([i[1] for i in location_prob]))
		# check that all is well
		if sum([i[1] for i in location_prob]) != 1:
			raise Exception('locaiton probs do not sum to 1')
		axins = ax.inset_axes(
				[node.absoluteTime-xSpan*0.1, node.y-ySpan*0.1, 
					xSpan*0.2, ySpan*0.2], 
				transform=ax.transData)
		axins.set_xticks([])
		axins.set_yticks([])
		axins.pie([i[1] for i in location_prob], 
			colors=[plot_config['colors'][i[0]] for i in location_prob],
			radius=0.5, frame=False,
			wedgeprops={"edgecolor":plot_config['colors']['base'],'linewidth': 1},
			normalize=False)
		ax.text(node.absoluteTime - xSpan*0.07, node.y-ySpan*0.07, 
			round(node.traits['posterior'], 2), 
			color=plot_config['colors']['base'])
	return(ax)


def plot_annotated_beast(ax, tree_file=None, plot_config=None, annotate_mrca_file=None): 
	myTree = import_tree(tree_file, plot_config)
	if annotate_mrca_file:
		# read in tip names to focus on
		focus_tip_names = \
			set(pd.read_csv(annotate_mrca_file, header=None)[0].str.replace('-', '_').str.split('|', expand=True)[1])
		# get the actual focus tips
		focus_tips = \
			[i for i in myTree.Objects if 
				i.branchType == 'leaf' and i.name.split('|')[1] in focus_tip_names]
		nodes_to_annotate = [myTree.commonAncestor(focus_tips), 
			myTree.commonAncestor(focus_tips).parent, 
			myTree.commonAncestor(focus_tips).parent.parent]
		print(nodes_to_annotate)
	else:
		nodes_to_annotate = [myTree.Objects[0]]
	xSpan, xLims, ySpan, yLims = \
		get_lims(myTree)
	ax = myTree.plotTree(ax, 
		x_attr=lambda k: k.absoluteTime,
		width=1.25,
		colour=lambda k: plot_config['colors'][k.traits['location']] if 
			k.traits['location'] in plot_config['labels'].keys() else 
				plot_config['colors']['Other (USA)'] if 'USA' in k.traits['location'] 
					else  plot_config['colors']['Other (intl.)'])
	ax = myTree.plotPoints(ax, 
		x_attr=lambda k: k.absoluteTime, 
		size=15, 
		colour=lambda k: plot_config['colors'][k.traits['location']] if 
			k.traits['location'] in plot_config['labels'].keys() else 
				plot_config['colors']['Other (USA)'] if 'USA' in k.traits['location'] 
					else  plot_config['colors']['Other (intl.)'],
		outline_colour=plot_config['colors']['base'])
	ax.set_xlim(xLims[0]-xSpan*0.1, xLims[1]+xSpan*0.1)
	ax.set_ylim(yLims[0]-ySpan*0.1, yLims[1]+ySpan*0.1)
	ax = \
		annotate_nodes(ax, myTree, xSpan, ySpan, nodes_to_annotate,  plot_config)
	ax_height = (ax.get_position().y1 - ax.get_position().y0) * ax.figure.bbox.height
	font_size = 10+int(ax.figure.bbox.height*0.01)*1.5
	starting_y = (1.3*font_size/ax_height)*(len(plot_config['labels'].items()))
	for label_idx, (label_key, label) in enumerate(plot_config['labels'].items()):
		text = ax.text(0, starting_y-(1.3*font_size/ax_height)*label_idx, 
			label, color=plot_config['colors'][label_key],
			size=font_size, transform=ax.transAxes, 
			path_effects=plot_config['effects'])
	return(ax, nodes_to_annotate)


def plot_table(ax, table_file=None, plot_config=None):
	x_vals = [0, 0.4, 0.7]
	y_vals = [0.75, 0.69]
	table_df = pd.read_csv(table_file, sep='\t', header=None)
	font_size = 10+int(ax.figure.bbox.height*0.01)*0.5
	# table headers
	ax.text(x_vals[0], y_vals[0], 'POS', 
		size=font_size, 
		fontweight='bold', va='top')
	ax.text(x_vals[1], y_vals[0], 'REF', 
		size=font_size, 
		fontweight='bold', va='top')
	ax.text(x_vals[2], y_vals[0], 'ALT', 
		size=font_size, 
		fontweight='bold', va='top')
	# table data
	ax.text(x_vals[0], y_vals[1], '\n'.join([str(i) for i in table_df[0].tolist()]),
		size=font_size,
		va='top')
	ax.text(x_vals[1], y_vals[1], '\n'.join([str(i) for i in table_df[1].tolist()]),
		size=font_size,
		va='top')
	ax.text(x_vals[2], y_vals[1], '\n'.join([str(i) for i in table_df[2].tolist()]),
		size=font_size,
		va='top')
	[ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left', 'bottom']]
	ax.set_xticks([])
	ax.set_yticks([])
	return(ax)


def plot_count_seqs_time(ax, match_profile_file=None, 
			loc_col=1, date_col=2, 
			cases_col=6, n_locs=6, plot_config=None, focal_region='USA'):
	match_seqs = pd.read_csv(match_profile_file, sep='\t', header=None)
	# adds week column
	match_seqs[date_col] = pd.to_datetime(match_seqs[date_col])
	match_seqs['week'] = match_seqs[date_col].dt.isocalendar().week
	# tabulates number per week
	n_per_week = match_seqs.groupby([loc_col, 'week']).size().reset_index()
	# sums by geo and gets top N geos
	top_n_locs = \
		set(n_per_week[[loc_col, 0]].groupby([loc_col]).sum().reset_index().sort_values(by=0,ascending=False)[:n_locs][loc_col])
	# add new location column
	new_location = \
		match_seqs.apply(lambda k: k[loc_col] if k[loc_col] in 
			top_n_locs else 'Other (USA)' if 'USA' in k[loc_col] else 'Other (intl.)', 
			axis=1).copy(deep=True)
	match_seqs['new_location'] = new_location
	n_per_week = \
		match_seqs.groupby(['week', 'new_location']).size().reset_index()
	n_per_week = n_per_week.pivot(index='week', columns='new_location', values=0)
	n_per_week.index = pd.to_datetime(n_per_week.index.astype(str)+
                           '2020-1',format='%V%G-%u')
	col_order = n_per_week.sum(axis=0).sort_values(ascending=False).index.tolist()
	col_order.remove('Other (USA)')
	col_order.remove('Other (intl.)')
	col_order.extend(['Other (USA)', 'Other (intl.)'])
	n_per_week = n_per_week[col_order]
	cols = [plot_config['colors'][i] for i in col_order]
	cm = LinearSegmentedColormap.from_list('states', cols, N=len(cols))
	n_per_week.plot.area(ax=ax, colormap=cm, alpha=0.95)
	ax.set_ylabel('Sequences/Wk')
	ax.set_xlabel('Date (2020)')
	font_size = 10+int(ax.figure.bbox.height*0.01)*1.35
	ax_height = (ax.get_position().y1 - ax.get_position().y0) * ax.figure.bbox.height
	for col_idx, col in enumerate(col_order[0:4]):
		ax.text(0.53, 0.85-(1.35*font_size/ax_height)*col_idx, 
			col, color=plot_config['colors'][col], 
			transform=ax.transAxes,
			size=font_size, 
			path_effects=plot_config['effects'])
	for col_idx, col in enumerate(col_order[4: ]):
		ax.text(0.76, 0.85-(1.35*font_size/ax_height)*col_idx, 
			col, color=plot_config['colors'][col], 
			transform=ax.transAxes,
			size=font_size, 
			path_effects=plot_config['effects'])
	ax.get_legend().remove()
	x_ticks = [
		 datetime.strptime('2020-03-01', '%Y-%m-%d'),
		 datetime.strptime('2020-03-15', '%Y-%m-%d'),
		 datetime.strptime('2020-04-01', '%Y-%m-%d'),
		 datetime.strptime('2020-04-15', '%Y-%m-%d'),
		 datetime.strptime('2020-05-01', '%Y-%m-%d'),
		 datetime.strptime('2020-05-15', '%Y-%m-%d'),
		 datetime.strptime('2020-06-01', '%Y-%m-%d'),
		 datetime.strptime('2020-06-15', '%Y-%m-%d'),
		 datetime.strptime('2020-07-01', '%Y-%m-%d'),
		 datetime.strptime('2020-07-15', '%Y-%m-%d'),
		 datetime.strptime('2020-08-01', '%Y-%m-%d')]
	x_labels = \
        [i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_labels, rotation=0, ha='center')
	ax.grid(axis='x', ls='-', color='#d8dee9', zorder=0)
	ax.set_axisbelow(True)
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
	parser.add_argument('--treeFocusTips',
		help='list of sequences to focus on')
	parser.add_argument('--mutationalProfile', 
						help='file with the mutational profile of sequences in tree')
	parser.add_argument('--matchMutationalProfile',
		help='estimated number of cases')
	parser.add_argument('--matchLocCol',
		help='which column has country data', 
		type=int, 
		default=0)
	parser.add_argument('--matchDateCol',
		help='which column has week numbers',
		type=int,
		default=0)
	parser.add_argument('--config',
		help='json file with config information (labels, colors, etc.)')
	args = parser.parse_args()
	plot_style()
	#args.tree = 'data/19B_subclade/19B_location_mcc_median.tre'
	#args.mutationalProfile = 'data/19B_subclade/19B_subclade_EHC_focal_snps_shared.tsv'
	#args.config = 'config/annotated_beast_tree_new2.json'
	#args.matchMutationalProfile = 'data/l84s/l84s_19B_subclade_match_focal_snps.tsv'
	args.matchLocCol = 8
	args.matchDateCol = 3
	#args.treeFocusTips = 'data/weighted_downsampling_old/ga_focused_aligned_masked_weighted_cluster_0.tsv'
	# reads in and processes config file
	plot_config = json.load(open(args.config, 'r'))
	plot_config['effects'] = eval(plot_config['effects'])
	plot_config['max_date'] = float(plot_config['max_date'])

		# todo generate programmaticaly
		# or at least read from config file
	x_ticks = [
		datetime.strptime('2020-01-01', '%Y-%m-%d'),
		datetime.strptime('2020-01-15', '%Y-%m-%d'),
		 datetime.strptime('2020-02-01', '%Y-%m-%d'),
		 datetime.strptime('2020-02-15', '%Y-%m-%d'),
		 datetime.strptime('2020-03-01', '%Y-%m-%d'),
		 datetime.strptime('2020-03-15', '%Y-%m-%d'),
		 datetime.strptime('2020-04-01', '%Y-%m-%d')]
	x_labels = \
        [i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
	x_ticks = [numeric_date(i) for i in x_ticks]


	fig = plt.figure(figsize=(6.4*1.59,4.8*2.09))
	gs=GridSpec(8,4, figure=fig)
	ax0 = fig.add_subplot(gs[:5,:3])
	ax1 = fig.add_subplot(gs[:5,3])
	ax2 = fig.add_subplot(gs[5:,:])
	axs = [ax0, ax1, ax2]
	axs[0], annotated_nodes = plot_annotated_beast(axs[0], tree_file=args.tree, 
		plot_config=plot_config, annotate_mrca_file=args.treeFocusTips)
	axs[0].grid(axis='x', ls='-', color='#d8dee9')
	axs[0].set_xticks(x_ticks)
	axs[0].set_xticklabels(x_labels)
	axs[0].set_xlabel('Date (2020)')
	[axs[0].spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
	axs[0].set_yticks([])
	axs[1] = plot_table(axs[1], table_file=args.mutationalProfile, plot_config=plot_config)
	axs[2] = \
		plot_count_seqs_time(axs[2], match_profile_file=args.matchMutationalProfile, 
			loc_col=args.matchLocCol, date_col=args.matchDateCol, 
			n_locs=6,  plot_config=plot_config)
	fig.tight_layout(h_pad=5.0)

	y_poses = [1.0, 1.0, 1.1]
	for ax_idx, ax in enumerate(axs):
	    ax.text(-0.09, y_poses[ax_idx], string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
	            size=20, weight='bold', va="top")
	    ax.set_axisbelow(True)

	fig.tight_layout()
	fig.savefig(f'{plot_config["out_name"]}.pdf')
	plt.close()



	# saves annotated node information
	with open(f'{plot_config["out_name"]}_annotated_nodes.txt', 'w') as fp:
		for node in annotated_nodes:
			fp.write(str(node))
			fp.write('\n')
			for key, value in vars(node).items():
				fp.write(str(key))
				fp.write(':\n')
				fp.write(str(value))
				fp.write('\n')
			fp.write('\n')


if __name__ == "__main__":
    run()

