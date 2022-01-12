import argparse
import pandas as pd 
from Bio import Phylo


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--tree', 
	    help='newick file with tree on which to cluster sequences from')
	parser.add_argument('--treeNameSep',
	    default='|',
	    help='character to seperate newick tree names on to match input seq names')
	parser.add_argument('--treeNameField',
	    default=1,
	    type=int,
	    help='which field in the character seperated tree names to take to match input seq names')
	parser.add_argument('--tipNames',
	    help='newline delimited text file of sequence names to cluster')
	parser.add_argument('--tipNameSep',
	    help='character to seperate sequence names',
	    default='|')
	parser.add_argument('--tipNameField',
	    help='field to take in seperated sequence names',
	    type=int,
	    default=1)
	parser.add_argument('--nGenerations', 
		type=int, 
		help='number of generations back to go',
		default=0)
	args = parser.parse_args()
	# reads in tree
	tree = Phylo.read(args.tree, 'newick')
	# generate tip dict
	tip_dict = {i.name.split(args.treeNameSep)[args.treeNameField]: i for i in tree.get_terminals()}
	# reads in seq names
	tip_names = set(pd.read_csv(args.tipNames, 
		header=None)[0].str.split(args.tipNameSep, expand=True)[args.tipNameField])
	# get MRCA of these sequences
	mrca = tree.common_ancestor([tip_dict[i] for i in tip_names])
	for generation in range(args.nGenerations):
		# check if we're already at the root
		if mrca != tree.get_nonterminals()[0]:
			new_mrca = tree.get_path(mrca)[-2]
			mrca = new_mrca

	descendant_tips = [i.name for i in mrca.get_terminals()]
	pd.DataFrame(descendant_tips).to_csv(args.tipNames.split('.')[0]+'_family.tsv', header=None, index=None) 


if __name__ == "__main__":
    run()


