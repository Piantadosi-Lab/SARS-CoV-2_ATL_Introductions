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
import datetime
import pandas as pd
from utils import plot_style, datetime_from_numeric
import seaborn as sns

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
		print(datetime_from_numeric(plot_config['max_date'] - node.traits["height"]))
		xmin = plot_config['max_date'] - node.traits["height_95%_HPD"][1]
		xmax = plot_config['max_date'] - node.traits["height_95%_HPD"][0]
		print(datetime_from_numeric(xmin))
		print(datetime_from_numeric(xmax))
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
		axins.text(-0.8, -0.4, 
			round(node.traits['posterior'], 1), 
			color=plot_config['colors']['base'], size=10, va='center', ha='center')
		print(node.traits['posterior'])
		print(round(node.traits['posterior'], 1))
		#ax.text(node.absoluteTime - xSpan*0.07, node.y, 
		#	round(node.traits['posterior'], 1), 
		#	color=plot_config['colors']['base'], size=10, va='center')
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
		nodes_to_annotate = [myTree.commonAncestor(focus_tips).parent, 
			myTree.commonAncestor(focus_tips).parent.parent,
			myTree.commonAncestor(focus_tips).parent.parent.parent]
			#myTree.commonAncestor(focus_tips).parent.parent]
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
	parser.add_argument('--config',
		help='json file with config information (labels, colors, etc.)')
	parser.add_argument('--nImportations')
	args = parser.parse_args()
	plot_style()
	#args.tree = 'data/19B_subclade/19B_location_mcc_median.tre'
	#args.mutationalProfile = 'data/19B_subclade/19B_subclade_EHC_focal_snps_shared.tsv'
	#args.config = 'config/annotated_beast_tree_new2.json'
	#args.matchMutationalProfile = 'data/l84s/l84s_19B_subclade_match_focal_snps.tsv'
	#args.matchLocCol = 8
	#args.matchDateCol = 3
	#args.treeFocusTips = 'data/weighted_downsampling_old/ga_focused_aligned_masked_weighted_cluster_0.tsv'
	# reads in and processes config file
	plot_config = json.load(open(args.config, 'r'))
	plot_config['effects'] = eval(plot_config['effects'])
	plot_config['max_date'] = float(plot_config['max_date'])

		# todo generate programmaticaly
		# or at least read from config file
	x_ticks = [
		datetime.datetime.strptime('2020-01-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-01-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-02-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-02-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-03-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-03-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-04-01', '%Y-%m-%d')]
	x_labels = \
        [i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
	x_ticks = [numeric_date(i) for i in x_ticks]

	from collections import Counter
	n_importations = Counter(pd.read_csv(args.nImportations, header=None)[0].values)
	total = sum(n_importations.values())
	print(total)
	n_importations = {key: value/total for key, value in n_importations.items()}
	print('HERE')
	print(n_importations)
	fig, axs = plt.subplots(1, 2, figsize=(6.4*1.25, 4.8), 
		constrained_layout=True, gridspec_kw={'width_ratios': [0.7, 0.3]})
	ax = axs[0]
	ax, annotated_nodes = plot_annotated_beast(ax, tree_file=args.tree, 
		plot_config=plot_config, annotate_mrca_file=args.treeFocusTips)
	ax.set_axisbelow(True)
	ax.grid(axis='x', ls='-', color='#d8dee9')
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_labels)
	ax.set_xlabel('Date (2020)')
	[ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
	ax.set_yticks([])
	axs[1].bar(n_importations.keys(), n_importations.values(), 
		facecolor=plot_config['colors']['GeorgiaUSA'], edgecolor='#333333', width=1.0)
	axs[1].set_xlabel('Number of introductions')
	axs[1].set_ylabel('Proportion')
	axs[1].set_aspect(1./axs[1].get_data_ratio())
	_ = [axs[1].spines[loc].set_visible(False) for loc in ['top', 'right']]
	x_poss = [0, -0.3]
	y_poss = [1.0, 1.49]
	for ax_idx, ax in enumerate(axs):
		ax.text(x_poss[ax_idx], y_poss[ax_idx], string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
			size=20, weight='bold', va="top")

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

