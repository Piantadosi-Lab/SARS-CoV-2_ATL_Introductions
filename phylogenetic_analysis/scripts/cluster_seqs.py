import argparse
import pandas as pd 
from Bio import Phylo
import itertools
import numpy as np
import scipy.sparse



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
    parser.add_argument('--seqNames',
        help='newline delimited text file of sequence names to cluster')
    parser.add_argument('--threshold',
        help='distance threshold to split clusters',
        type=float,
        default=0.3)
    args = parser.parse_args()
    #args.seqNames = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv'
    #args.tree = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick'
    tree = Phylo.read(open(args.tree, 'r'), 'newick')
    # biopython does not read in branch lengths correctly, weird hack -- revisit todo
    #for node in tree.find_clades():
    #    node.branch_length = float(node.branch_length.lstrip('_').replace('_', '-'))

    tip_dict = {i.name.split(args.treeNameSep)[args.treeNameField]: i for i in tree.get_terminals()}
    if args.seqNames:
        get_tip_names = set(pd.read_csv(args.seqNames, sep='\t', header=None)[1])
        get_tips = [tip_dict[i] for i in get_tip_names if i in tip_dict.keys()]
        print([i for i in get_tip_names if i not in tip_dict.keys()])
        print(f'{len(get_tips)} tips in --seqNames file found in tree')

    else:
        get_tips = list(tip_dict.values())

    # use permutation so matrix is symmetrical
    pairs = itertools.permutations(get_tips,2)
    dists = \
        [[pair[0].name, pair[1].name, tree.distance(pair[0], pair[1])] for pair in pairs]
    dists = pd.DataFrame(dists)
    dists[3] = dists[2].apply(lambda x: x < args.threshold).astype(int)
    # generate matrix
    mat_df = dists[[0,1,3]].pivot(index=0, columns=1)
    mat = np.array(mat_df)
    cc = scipy.sparse.csgraph.connected_components(mat)
    cc_assigned = pd.DataFrame(zip(mat_df.index, cc[1]))

    print(f'there are {len(cc_assigned[1].unique())} clusters using a threshold of {args.threshold}')
    for cc_idx, cc_id in enumerate(np.unique(cc[1])):
        print(f'there are {np.unique(cc[1], return_counts=True)[1][cc_idx]} sequences in cluster {cc_id}')

    cc_assigned.to_csv(args.tree.split('.')[0]+f'_clusters_{args.threshold}.tsv', sep='\t', header=None, index=None)


if __name__ == "__main__":
    run()