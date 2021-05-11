import pickle
import argparse
import pandas as pd
import numpy as np

def get_importations(tree_idx, tree, translate_dict, focal_region):
	parse_loc = lambda k: k.comment.split('=')[1].replace('"', '')
	# get all tips from focal region
	focal_tips = [i for i in tree.get_terminals() if 
		parse_loc(i) == focal_region]
	# get the path from root to tip for all focal tips
	# root is not automatically returned so we add it in
	focal_paths = [[tree.root, *tree.get_path(item)] for item in focal_tips]
	# get the location for each node in the path
	focal_paths_locs = \
		[[parse_loc(item) for item in path] 
			for path in focal_paths]
	# index of node which precede the earliest node assigned to the focal region
	# minimum value is 0, root of the tree
	importation_node_indices = \
		[min([i for i, item in enumerate(path) if item==focal_region])-1
		for path in focal_paths_locs]
	importation_node_indices = [i if i != -1 else np.nan for i in importation_node_indices]
	importation_nodes = [focal_paths[item][i] if not np.isnan(i) else tree.root for item, i in enumerate(importation_node_indices)]
	# reduces to unique nodes with dictionary hack
	importation_nodes = list({i.name: i for i in importation_nodes}.values())
	# for each importation node, get the earliest descendant assigned to focal_region
	# if root, then imported node == importation node
	imported_nodes = [sorted([item for item in i.clades if parse_loc(item)==focal_region], 
			key=lambda k: k.branch_length)[0] if i != tree.root else i for i in importation_nodes]
	importation_node_source = \
		[i.comment.split('=')[1].replace('"', '') for 
			i in importation_nodes]
	importation_node_height = \
		[tree.distance(i)+imported_nodes[idx].branch_length/2 for 
			idx, i in enumerate(importation_nodes)]
	importation_node_descendants = \
		[[item.name for item in i.get_terminals() if parse_loc(item)==focal_region]
			for i in importation_nodes]
	importation_node_descendants = \
		[[translate_dict[str(item)] for item in i] for i 
			in importation_node_descendants]
	importation_node_n_descendants = [len(i) for i in importation_node_descendants]
	importation_nodes_df = \
		pd.DataFrame([importation_node_source, importation_node_height, 
			importation_node_n_descendants, importation_node_descendants]).transpose()
	importation_nodes_df['tree_idx'] = tree_idx
	n_importations = len(importation_nodes)
	return(n_importations, importation_nodes_df)


# mostly from arviz
# but works directly with a list
# and also returns median
# https://github.com/arviz-devs/arviz/blob/master/arviz/stats/stats.py
def hpd(dat, qs):
    width = qs[2] - qs[0]
    # sorts from smallest to largest
    dat = sorted(dat)
    # length of data
    n = len(dat)
    # number of values we are keeping
    # thus, the index of the begining of 
    # the HPD interval must be <= this value
    # this gives us the tail of the distribution
    interval_idx_inc = int(np.floor(width * n))
    # number of values we are excluding
    # thus, possible number of HPD intervals
    # this gives us the head of the distribution
    n_intervals = n - interval_idx_inc
    # the width of each possible interval
    # for each possible head and tail value, 
    # what is the difference between them
    interval_width = [a_i - b_i for a_i, b_i in 
                      zip(dat[interval_idx_inc:], 
                          dat[:n_intervals])]
    # find the shortest interval
    min_idx = interval_width.index(min(interval_width))
    hpd_interval = (dat[min_idx], dat[min_idx+interval_idx_inc])
    dat_hpd = [item for item in dat if (item >= hpd_interval[0]) & (item <= hpd_interval[1])]
    dat_mid = np.quantile(dat_hpd, qs[1])
    return((hpd_interval[0], dat_mid, hpd_interval[1]))


def run():
	# parses beast nexus trees into biopython tree objects
	# having trouble with BioPython nexus parser
	# todo revisit whether this is necessary
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--trees', 
		help='path to set of parsed beast trees, assumed pickle format')
	parser.add_argument('--focalRegion',
		help='which region to estimate importations into')
	args = parser.parse_args()
	#args.trees = 'data/19B_subclade/19B_sublcade_location_tree_with_trait.pkl'
	#args.focalRegion = 'GeorgiaUSA'
	translate_dict, trees = pickle.load(open(args.trees, 'rb'))
	print(len(trees))
	importations_df = pd.DataFrame()
	n_importations = []
	for tree_idx, tree in enumerate(trees):
		tree_importations = get_importations(tree_idx, tree, translate_dict, args.focalRegion)
		n_importations.append(tree_importations[0])
		importations_df = pd.concat([importations_df, tree_importations[1]])
	print(importations_df)
	print(hpd(np.array(n_importations), [0.025, 0.5, 0.975]))
	pd.DataFrame(n_importations).to_csv('.'.join(args.trees.split('.')[0:-1]) + 
		'_importations_count.tsv', sep='\t', header=None, index=None)
	importations_df.to_csv('.'.join(args.trees.split('.')[0:-1]) + 
		'_importations_dat.tsv', sep='\t', header=None, index=None)




if __name__ == '__main__':
    run()