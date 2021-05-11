#!/usr/bin/env python3
import subprocess
import shlex
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import treetime as tt
from treetime import TreeTime
from treetime import wrappers
from treetime.utils import numeric_date
import pandas as pd
import argparse
import ast
import numpy as np
import os
from collections import defaultdict
import io
import sys


def label_nodes(tre):
	node_counter=0
	for item in tre.get_nonterminals():
		if item.name is None:
			item.name = 'NODE_'+format(node_counter, '07d')
			node_counter+=1
	return(tre)


def infer_timetree(tre, aln_path, date_dict, resolve, outgroup, save_filtered_aln, out_path):
	# dates are for ClockTree class, all other ClockTree parameters are default
	# tree, aln, gtr arguments are passed to TreeAnc clas
	tre = TreeTime(dates=date_dict, tree=tre, aln=aln_path, gtr='JC69')
	# from NextStrain Augur refine.py
	# treetime clock filter will mark, but not remove bad tips
	tre.clock_filter(reroot=outgroup, n_iqd=4, plot=False) 
	# remove them explicitly
	leaves = [x for x in tre.tree.get_terminals()]
	bad_seqs = [item.name for item in leaves if item.bad_branch==True]
	pd.DataFrame(bad_seqs).to_csv(f'{out_path}_bad_seqs.tsv', sep='\t', header=None, index=None)
	print(f'{len(bad_seqs)} tips failed clock filter')
	seqs = list(SeqIO.parse(aln_path, 'fasta'))
	seqs = [item for item in seqs if item.id not in bad_seqs]
	if save_filtered_aln:
		time_filtered_seqs_path = out_path+'_time.fasta'
		with open(time_filtered_seqs_path, 'w') as out:
			SeqIO.write(seqs, out, 'fasta')
	for n in leaves:
		if n.bad_branch:
			tre.tree.prune(n)
			print('pruning leaf ', n.name)
	# fix treetime set-up for new tree topology
	tre.prepare_tree()
	print(f'resolve polytomies: {resolve}')
	tre.run(root=outgroup, max_iter=2, tc='skyline', 
		fixed_clock_rate=0.001, time_marginal=True, 
		vary_rate=0.0005, resolve_polytomies=resolve)
	# save data for root to tip regression
	rtt_dat = [[n.numdate, n._v] for n in tre.tree.get_terminals()]
	pd.DataFrame(rtt_dat).to_csv(f'{out_path}_rtt.csv', header=None, index=None)
	times = pd.DataFrame({'name': [item.name for item in tre.tree.find_clades()],
						  'date': [item.numdate for item in tre.tree.find_clades()],
						  'lower': [list(tre.get_max_posterior_region(item, 0.9))[0] 
						  			for item in tre.tree.find_clades()],
						  'upper': [list(tre.get_max_posterior_region(item, 0.9))[1] 
						  			for item in tre.tree.find_clades()]},
						  index = range(0, len([item for item in tre.tree.find_clades()])))
	times.to_csv(f'{out_path}_refined_node_times.csv')
	for n in tre.tree.find_clades():
		n.branch_length = n.mutation_length
	# saving as XML to avoid recurssion error
	with open(f'{out_path}_refined.newick', 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'newick')

	tre.branch_length_to_years()
	with open(f'{out_path}_refined_time.newick', 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'newick')
	return(tre)


def infer_mugration(tre, aln_file, loc_dict, bias_correction, out_path):
	# associates tree/seq names with metadata names
	# assumes tree and seq names match
	# metadata has more rows that tips in tree
	# remove amp and merged is hack, fix later
	tre, letter_to_state, state_to_letter = \
		wrappers.reconstruct_discrete_traits(tre, loc_dict, 
			sampling_bias_correction=bias_correction)
	# saves gtr model
	# from treetiem mugration code
	with open(out_path+'_gtr.txt', 'w', encoding='utf-8') as ofile:
		ofile.write('Character to attribute mapping:\n')
		for state in letter_to_state.values():
			ofile.write('  %s: %s\n'%(state_to_letter[state], state))
		ofile.write('\n\n'+str(tre.gtr)+'\n')
		print("\nSaved inferred mugration model")
	region_letters = pd.DataFrame({'letter':list(letter_to_state.keys()), 
		'region': list(letter_to_state.values())})
	region_letters.to_csv(f'{out_path}_refined_regions.csv', index=False)
	node_states = pd.DataFrame({'name': [item.name for item in tre.tree.find_clades()],
		'state': [letter_to_state[item.cseq[0]] for item in tre.tree.find_clades()],
		'state_conf': [list(item.marginal_profile[0]) for item in tre.tree.find_clades()]})
	node_states.to_csv(f'{out_path}_refined_node_states.csv', index=False)
	return(tre, letter_to_state, state_to_letter)


def estimate_importations(tre, loc_dict, interest_regions, out_path):
	# list of regions and their identifying character
	region_dict = pd.read_csv(f'{out_path}_refined_regions.csv')
	region_dict = {i['letter']: i['region'] for idx, i in region_dict.iterrows()}
	# Dictionary of each node's state
	node_states = pd.read_csv(f'{out_path}_refined_node_states.csv')
	node_states_dict = \
		{item['name']: item['state'] for index, item in node_states.iterrows()}
	# Dictionary of each node's time
	node_times = \
		pd.read_csv(f'{out_path}_refined_node_times.csv')
	node_times_dict = \
		{item['name']: (item['lower'], item['date'], item['upper']) for index, item in node_times.iterrows()}
	# formats the unique importation nodes list
	formatted_unique_importation_nodes = []
	formatted_columns = ['destination region', 'source region', 'descendant tips', 
					     'source node', '# descendants', 'midroot time', 'lower midroot time', 'upper midroot time']
	# tabulates the number of introductions into each region
	n_introductions_regions = {}
	for region in interest_regions:
		# List of tree tips from specific region of interest
		region_tre_tips = \
			[item for item in tre.get_terminals() if node_states_dict[item.name] == region]
		# For each tip from the region of interest, get path from root to ip
		region_paths = [tre.get_path(item) for item in region_tre_tips]
		# state of each step in each path
		region_paths_states = [[node_states_dict[item.name] for item in path] for path in region_paths]
		# node which precede the earliest node assigned to ANY of the regions of interest
		importation_node_indices = \
			[min([i for i, item in enumerate(path) if item in interest_regions])-1 
			for path in region_paths_states]
		importation_nodes = \
			[region_paths[i][item] for i, item in enumerate(importation_node_indices)]
		importation_nodes_dict = \
			{region_paths[i][-1]: item for i, item in enumerate(importation_nodes)}
		importation_nodes = list(set(importation_nodes))
		n_introductions_regions[region] = len(importation_nodes)
		for importation_node in importation_nodes:
			node_time = node_times_dict[importation_node.name]
			node_descendants = [key.name for key, value in importation_nodes_dict.items() 
								if value == importation_node]
			# sorts by the ML time
			min_descendant_time = sorted([node_times_dict[item] for item in node_descendants], 
										  key=lambda k: k[1])[0]
			this_importation_node = \
					[region, 
					 node_states_dict[importation_node.name], 
					 node_descendants, 
					 importation_node.name,
					 len(node_descendants),
					 (min_descendant_time[1] - node_time[1])/2 + node_time[1],
					 (min_descendant_time[0] - node_time[0])/2 + node_time[0],
					 (min_descendant_time[2] - node_time[2])/2 + node_time[2]]
			# adds output row to larger list
			formatted_unique_importation_nodes.append(this_importation_node)
	# turns formatted list into data frame
	formatted_unique_importation_nodes = \
		pd.DataFrame(formatted_unique_importation_nodes, columns=formatted_columns)
		# saves formatted list
	formatted_unique_importation_nodes.to_csv(f'{out_path}_refined_importations.csv')
	len(formatted_unique_importation_nodes['source node'])
	n_introductions = len(set(formatted_unique_importation_nodes['source node']))
	print(f'Estimated {n_introductions} introductions into all regions of interest')
	for region in interest_regions:
		n_region_introductions = \
			formatted_unique_importation_nodes[
				formatted_unique_importation_nodes['destination region'] == region].shape[0]
		print(f'Estimated {n_region_introductions} introductions into {region}')
	return(formatted_unique_importation_nodes)


# save gtr from states
def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--bootstrapReps', 
						help='index limits for bootstrap replicates to analyze',
						type=int,
						nargs=2,
						default=[0, 1])
	parser.add_argument('--trees', 
		help='file which contains set of trees to analyze',
		default='results/subsampled_alignment_neighbors.fasta.ufboot')
	parser.add_argument('--treeNameSep',
	    default='|')
	parser.add_argument('--treeNameField',
	    default=1, type=int, 
	    help='which field in the tree name corresponds to the metadata name column')
	parser.add_argument('--metadata',
						help='metadata file',
						default=None)
	parser.add_argument('--metadataDelim',
						help='metadata file delimiter',
						default='\t')
	parser.add_argument('--metadataNameCol',
						help='metadata column with squence names',
						type=int,
						default=1)
	parser.add_argument('--metadataDateCol',
						help='metadata column with sampling dates',
						type=int,
						default=2)
	parser.add_argument('--metadataDateFmt',
						help='date format',
						default='%Y-%m-%d')
	parser.add_argument('--metadataLocCol',
						help='metadata column with sampling location',
						type=int,
						default=6)
	parser.add_argument('--sequences', 
						help='alignment', 
						default=None)
	parser.add_argument('--outgroup',
					help='Outgroup to use in tree',
					default='Best')
	parser.add_argument('--biasCorrection',
					help='sampling bias correction',
					default=1, type=float)
	parser.add_argument('--regions',
					help='regions to estimate introductions into',
					default=None,
					nargs='+')
	parser.add_argument('--resolve_polytomies', dest='resolve', action='store_true')
	parser.add_argument('--saveAln', dest='save_aln', action='store_true')
	parser.set_defaults(resolve=False)
	parser.set_defaults(save_aln=False)
	args = parser.parse_args()
	sys.setrecursionlimit(5000)
	print('args parsed')
	# read in trees
	with open(args.trees, 'r') as fp:
		bootstrap_tre_strings = fp.readlines()
	print('trees read')
	# set up tree directories
	tre_rep_dir = f'{args.trees}_tres'
	if not os.path.exists(tre_rep_dir):
		os.mkdir(tre_rep_dir)
	# read in and format matadata
	metadata = pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
	if args.metadataDateFmt != 'numeric':
		metadata[args.metadataDateCol] = \
			pd.to_datetime(metadata[args.metadataDateCol], format=args.metadataDateFmt)
		metadata['numeric_date'] = metadata[args.metadataDateCol].apply(numeric_date)
	else:
		metadata['numeric_diate'] = metadata[args.metadataDateCol]
	# associates tree/seq names with metadata names
	# assumes tree and seq names match
	# all bootstrap replicates have same tips
	# remove merged and amp bit of a hack for now, remove from alignment names
	# loop over bootstrap replicates
	for i in range(*args.bootstrapReps):
		out_path = f'{args.trees}_tres/{i}'
		out_file = out_path + f'/{i}'
		if not os.path.exists(out_path):
			os.mkdir(out_path)
		print(f'replicate {i}')
		bootstrap_tre_string = bootstrap_tre_strings[i]
		bootstrap_tre = \
			Phylo.read(io.StringIO(bootstrap_tre_string), 'newick')
		bootstrap_tre = label_nodes(bootstrap_tre)
		print(f'replicate {i} tree internal nodes labelled')
		name_dict = defaultdict(lambda: 'not_in_tree')
		name_dict.update({i.name.split(args.treeNameSep)[args.treeNameField]: 
			i.name for i in bootstrap_tre.get_terminals()})
		# associates tree names with states (dates, locations)
		date_dict = {name_dict[i[args.metadataNameCol]]: i['numeric_date']
			for idx, i in metadata.iterrows()}
		loc_dict = {name_dict[i[args.metadataNameCol]]: i[args.metadataLocCol]
			for idx, i in metadata.iterrows()}
		out_group = name_dict[args.outgroup]
		tt_tre = \
			infer_timetree(bootstrap_tre, args.sequences, 
				date_dict, args.resolve, out_group, args.save_aln, out_file)
		print(f'replicate {i} timetree inferred')
		tt_tre2, letter_to_state, state_to_letter = \
			infer_mugration(tt_tre.tree, args.sequences, loc_dict, args.biasCorrection,
				out_file)
		print(f'replicate {i} mugrations inferred')
		importations = \
			estimate_importations(tt_tre2.tree, loc_dict, args.regions, out_file) 
		#refined_tree = Phylo.read(f'{out_file}_refined_time.newick', 'newick')
		#importations = \
		#	estimate_importations(refined_tree, loc_dict, args.regions, out_file)
		print(f'replicate {i} importations estimated')


if __name__ == "__main__":
	run()


