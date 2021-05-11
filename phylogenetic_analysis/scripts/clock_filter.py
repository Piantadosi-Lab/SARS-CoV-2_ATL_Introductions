import argparse
import pandas as pd 
from treetime import TreeTime
from treetime.utils import numeric_date
from Bio import Phylo, SeqIO

def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--seqs', 
	    help='path to sequences')
	parser.add_argument('--tree', 
	    help='newick file with tree on which to cluster sequences from')
	parser.add_argument('--treeNameSep',
	    default='|',
	    help='character to seperate newick tree names on to match input seq names')
	parser.add_argument('--treeNameField',
	    default=1,
	    type=int,
	    help='which field in the character seperated tree names to take to match input seq names')
	parser.add_argument('--metadata', 
	    default=None, help='metadata file path, assumes no header')
	parser.add_argument('--metadataDelim',
	    default='\t', 
	    help='metadata delimiter')
	parser.add_argument('--metadataIDCol',
	    default=1,
	    type=int,
	    help='which column in the metadata file contains sequence names (corresponds to seqNameField)')
	parser.add_argument('--metadataDateCol',
	    default=2,
	    type=int,
	    help='which column in the metadata file contains sequece dates')
	parser.add_argument('--metadataDateFmt',
	    default='%Y-%m-%d',
	    help='format of the date column')
	parser.add_argument('--nIQD',
	    default=4,
	    type=int)
	args = parser.parse_args()
	out_path = args.tree.replace('.treefile', '')
	# creates tree tip name dict
	tre = Phylo.read(args.tree, 'newick')
	tip_name_dict = {i.name.split(args.treeNameSep)[args.treeNameField]: 
		i.name for i in tre.get_terminals()}
	# read in metadata
	metadata = pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
	# set up date dictionary
	if args.metadataDateFmt != 'numeric':
		metadata[args.metadataDateCol] = \
			pd.to_datetime(metadata[args.metadataDateCol], format=args.metadataDateFmt)
		metadata['numeric_date'] = metadata[args.metadataDateCol].apply(numeric_date)
	else:
		metadata['numeric_date'] = metadata[args.metadataDateCol]
	date_dict = \
		{tip_name_dict[i[args.metadataIDCol]]: i['numeric_date'] for idx, i in metadata.iterrows()}
	tre = \
		TreeTime(dates=date_dict, tree=args.tree, aln=args.seqs, gtr='JC69')
	tre.clock_filter(reroot="best", n_iqd=4, plot=False) 
	# actually removes bad tips
	bad_tips = [item.name for item in tre.tree.get_terminals() if item.bad_branch==True]
	for n in tre.tree.get_terminals():
		if n.bad_branch:
			tre.tree.prune(n)
			print('pruning leaf ', n.name)
	# save report with bad tips
	pd.DataFrame(bad_tips).to_csv(f'{out_path}_bad_tips.tsv', sep='\t', header=None, index=None)
	# save filtered rerooted tree
	with open(f'{out_path}_clockfilter.newick', 'w') as out_file:
		Phylo.write(tre.tree, out_file, 'newick')
	# save filtered alignment
	tre_names = set([i.name for i in tre.tree.get_terminals()])
	seqs = list(SeqIO.parse(args.seqs, 'fasta'))
	tre_seqs = \
		[item for item in seqs if item.description.replace(' ', '_') in tre_names]
	with open(f'{out_path}_clockfilter.fasta', 'w') as out:
		SeqIO.write(tre_seqs, out, 'fasta')
	# save RTT info
	rtt_dat = [[n.raw_date_constraint, tre.tree.distance(n)] for n in tre.tree.get_terminals()]
	pd.DataFrame(rtt_dat).to_csv(f'{out_path}_rtt.csv', header=None, index=None)


if __name__ == "__main__":
	run()


