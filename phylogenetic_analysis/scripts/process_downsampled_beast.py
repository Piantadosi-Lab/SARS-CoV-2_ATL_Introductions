import baltic as bt
import pandas as pd
import argparse
import glob


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--trees')
	parser.add_argument('--treeNameSep', default='|')
	parser.add_argument('--treeNameField', default=1, type=int)
	parser.add_argument('--metadata')
	parser.add_argument('--metadataNameCol', type=int, default=1)
	parser.add_argument('--metadataDelim', default='\t')
	args = parser.parse_args()
	#args.trees = 'data/19B_subclade/19B_subclade_sampled_5_rep_*_mcc.tre'
	#args.metadata = 'data/19B_subclade/test.tsv'
	metadata = pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
	get_tips = set(metadata[args.metadataNameCol].values)
	tree_paths = glob.glob(args.trees)
	out_dat = []
	for tree_path_idx, tree_path in enumerate(tree_paths):	
		tree =  bt.loadNexus(tree_path, absoluteTime=False)
		tree_get_tips = get_tips & \
			set([i.name.split(args.treeNameSep)[args.treeNameField] for 
				i in tree.getExternal()])
		mrca = tree.commonAncestor(tree.getExternal(lambda k: 
			k.name.split(args.treeNameSep)[args.treeNameField] in get_tips))
		out_dat.extend(list(zip([tree_path_idx]*len(mrca.traits['location.set']), 
			mrca.traits['location.set'], 
			mrca.traits['location.set.prob'])))


	out_dat = pd.DataFrame(out_dat)
	out_dat.to_csv('.'.join(args.trees.split('.')[:-1]) + '_parsed.tsv', 
		sep='\t', header=None, index=None)





if __name__ == "__main__":
    run()



