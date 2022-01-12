from Bio import Phylo
import argparse
import pickle
import math
import io

def label_nodes(tre):
	node_counter=0
	for item in tre.get_nonterminals():
		if item.name is None:
			item.name = 'NODE_'+format(node_counter, '07d')
			node_counter+=1
	return(tre)


def run():
	# parses beast nexus trees into biopython tree objects
	# having trouble with BioPython nexus parser
	# todo revisit whether this is necessary
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--trees', 
		help='path to set of beast trees, assumed nexus format')
	parser.add_argument('--burnin', 
		help='what proportion of trees to discard as burnin',
		default=10,
		type=int)
	args = parser.parse_args()
	#args.trees = 'data/19B_subclade/19B_sublcade_location_tree_with_trait.trees'
	# reads in trees, very slow for large #s of trees!!
	# resaves as pickle to speed things ups
	with open(args.trees, 'r') as fp:
		trees = fp.readlines()

	translate_range = [[idx for idx,i in enumerate(trees) if 'Begin trees;' in i][0]]
	translate_range.append([idx for idx, i in enumerate(trees) if ';' in i and idx > translate_range[0]][0])
	translate_dict = \
		[i.replace(',\n', '').lstrip().split(' ') for i in trees[translate_range[0]+2:translate_range[1]]]
	translate_dict = {i[0]: i[1] for i in translate_dict}
	newick_trees = ['(' + '('.join(i.split('(')[1:]) for i in trees[translate_range[1]+1:] if i[:3] != 'End']
	discard = int(math.ceil(len(newick_trees)*args.burnin/100))
	print(f'discarding {args.burnin}%, {discard} trees as burnin')
	newick_trees = newick_trees[discard:]
	# slow with many trees!
	biopython_trees = [Phylo.read(io.StringIO(i), 'newick') for i in newick_trees]
	biopython_trees = [label_nodes(i) for i in biopython_trees]
	parsed_trees = (translate_dict, biopython_trees)
	print(len(biopython_trees))
	pickle.dump(parsed_trees, open('.'.join(args.trees.split('.')[0:-1])+'.pkl', 'wb'))



if __name__ == '__main__':
    run()