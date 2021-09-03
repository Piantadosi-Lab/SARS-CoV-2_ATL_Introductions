import pandas as pd
import glob
import argparse





def run()
parser = argparse.ArgumentParser()
parser.add_argument('--inDir')
parser.add_argument('--metadata')
args = parser.parse_args()
args.inDir = 'data/traveler_analysis'
args.metadata = 'data/ga_travel_seqs.tsv'

metadata = pd.read_csv(args.metadata,  sep='\t', header=None)

all_output = ''
for row_idx, row in metadata.iterrows():
	# travel region
	print(row)
	n_travel_seqs = \
		pd.read_csv(f'data/travel_analysis/{row[1]}_{row[2]}_n_genotyped.tsv', 
			sep='\t', header=None).values[0][2]
	travel_data = pd.read_csv(f'data/travel_analysis/{row[1]}_{row[2]}_dists.tsv', sep='\t', header=None)
	row_seq = travel_data[travel_data[0] == row[1]]
	travel_data['dist_to_traveler'] = travel_data[3] - row_seq[3].values[0]
	min_travel_snps = min(travel_data[travel_data[0] != row[1]]['dist_to_traveler'].abs()) if travel_data[travel_data[0] != row[1]].shape[0] > 0 else np.nan
	if min_travel_snps != 0 and ~np.isnan(min_travel_snps):
		ancestral = True if travel_data[travel_data['dist_to_traveler'] == -1 * min_travel_snps].shape[0] > 0 else False
		descendant = True if travel_data[travel_data['dist_to_traveler'] == min_travel_snps].shape[0] > 0 else False
		if (ancestral == True and descendant == False):
			min_travel_snps_label = " (ancestral)"
			min_travel_snps_count = \
				travel_data[(travel_data[0] != row[1]) & 
				(travel_data['dist_to_traveler'].abs() == min_travel_snps)].shape[0]
		elif (ancestral == False and descendant == True):
			min_travel_snps_label = " (descendant)"
			min_travel_snps_count = \
				travel_data[(travel_data[0] != row[1]) & 
				(travel_data['dist_to_traveler'].abs() == min_travel_snps)].shape[0]
		elif (ancestral == True and descendant == True):
			min_travel_snps_label = " (ancestral, descendant)"
			min_travel_snps_count = \
				f'{travel_data[(travel_data[0] != row[1]) & (-1 * travel_data["dist_to_traveler"] == min_travel_snps)].shape[0]} (ancestral); '
			min_travel_snps_count += \
				f'{travel_data[(travel_data[0] != row[1]) & (travel_data["dist_to_traveler"] == min_travel_snps)].shape[0]} (descendant)'
	elif ~np.isnan(min_travel_snps):
		min_travel_snps_label = ''
		min_travel_snps_count = \
			travel_data[(travel_data[0] != row[1]) & (travel_data["dist_to_traveler"] == min_travel_snps)].shape[0]
	else:
		min_travel_snps_label = ''
		min_travel_snps_count = np.nan
	travel_output = f'{min_travel_snps}{min_travel_snps_label}	{min_travel_snps_count}	{n_travel_seqs}'
	# georgia
	n_ga_seqs = \
		pd.read_csv(f'data/travel_analysis/{row[1]}_GeorgiaUSA_n_genotyped.tsv', 
			sep='\t', header=None).values[0][2]
	ga_data = pd.read_csv(f'data/travel_analysis/{row[1]}_GeorgiaUSA_dists.tsv', sep='\t', header=None)
	row_seq = ga_data[ga_data[0] == row[1]]
	ga_data['dist_to_traveler'] = ga_data[3] - row_seq[3].values[0]
	min_ga_snps = min(ga_data[ga_data[0] != row[1]]['dist_to_traveler'].abs()) if ga_data[ga_data[0] != row[1]].shape[0] > 0 else np.nan
	if min_ga_snps != 0 and ~np.isnan(min_ga_snps):
		ancestral = True if ga_data[ga_data['dist_to_traveler'] == -1 * min_ga_snps].shape[0] > 0 else False
		descendant = True if ga_data[ga_data['dist_to_traveler'] == min_ga_snps].shape[0] > 0 else False
		if (ancestral == True and descendant == False):
			min_ga_snps_label = " (ancestral)"
			min_ga_snps_count = \
				ga_data[(ga_data[0] != row[1]) & 
				(ga_data['dist_to_traveler'].abs() == min_ga_snps)].shape[0]
		elif (ancestral == False and descendant == True):
			min_ga_snps_label = " (descendant)"
			min_ga_snps_count = \
				ga_data[(ga_data[0] != row[1]) & 
				(ga_data['dist_to_traveler'].abs() == min_ga_snps)].shape[0]
		elif (ancestral == True and descendant == True):
			min_ga_snps_label = " (ancestral, descendant)"
			min_ga_snps_count = \
				f'{ga_data[(ga_data[0] != row[1]) & (-1 * ga_data["dist_to_traveler"] == min_ga_snps)].shape[0]} (ancestral); '
			min_ga_snps_count += \
				f'{ ga_data[(ga_data[0] != row[1]) & (ga_data["dist_to_traveler"] == min_ga_snps)].shape[0]} (descendant)'
	elif ~np.isnan(min_ga_snps):
		min_ga_snps_label = ''
		min_ga_snps_count = \
				ga_data[(ga_data[0] != row[1]) & (ga_data["dist_to_traveler"] == min_ga_snps)].shape[0]
	else:
		min_ga_snps_label = ''
		min_ga_snps_count = np.nan
	ga_output = f'{min_ga_snps}{min_ga_snps_label}	{min_ga_snps_count}	{n_ga_seqs}'
	# finally, global
	n_global_seqs = \
		pd.read_csv(f'data/travel_analysis/{row[1]}_None_n_genotyped.tsv', 
			sep='\t', header=None).values[0][2]
	# need to make out GA and travel region
	n_global_seqs = n_global_seqs - n_ga_seqs - n_travel_seqs
	global_data = pd.read_csv(f'data/travel_analysis/{row[1]}_None_dists.tsv', sep='\t', header=None)
	row_seq = global_data[global_data[0] == row[1]].copy()
	# now back out GA and travel region seqs
	global_data = global_data[(global_data[1] != 'GeorgiaUSA') & (global_data[1] != row[2])]
	global_data['dist_to_traveler'] = global_data[3] - row_seq[3].values[0]
	min_global_snps = min(global_data[global_data[0] != row[1]]['dist_to_traveler'].abs()) if global_data[global_data[0] != row[1]].shape[0] > 0 else np.nan
	if min_global_snps != 0 and ~np.isnan(min_global_snps):
		ancestral = True if global_data[global_data['dist_to_traveler'] == -1 * min_global_snps].shape[0] > 0 else False
		descendant = True if global_data[global_data['dist_to_traveler'] == min_global_snps].shape[0] > 0 else False
		if (ancestral == True and descendant == False):
			min_global_snps_label = " (ancestral)"
			min_global_snps_count = \
				global_data[(global_data[0] != row[1]) & 
				(global_data['dist_to_traveler'].abs() == min_global_snps)].shape[0]
		elif (ancestral == False and descendant == True):
			min_global_snps_label = " (descendant)"
			min_global_snps_count = \
				global_data[(global_data[0] != row[1]) & 
				(global_data['dist_to_traveler'].abs() == min_global_snps)].shape[0]
		elif (ancestral == True and descendant == True):
			min_global_snps_label = " (ancestral, descendant)"
			min_global_snps_count = \
				f'{global_data[(global_data[0] != row[1]) & (-1 * global_data["dist_to_traveler"] == min_global_snps)].shape[0]} (ancestral); '
			min_global_snps_count += \
				f'{global_data[(global_data[0] != row[1]) & (global_data["dist_to_traveler"] == min_global_snps)].shape[0]} (descendant)'
	elif ~np.isnan(min_global_snps):
		min_global_snps_label = ''
		min_global_snps_count = global_data[(global_data[0] != row[1]) & (global_data["dist_to_traveler"] == min_global_snps)].shape[0]
	else:
		min_global_snps_label = ''
		min_global_snps_count = np.nan
	min_global_snps_locs = ';'.join([i.replace('USA', ', USA') for i in np.unique(global_data[(global_data[0] != row[0]) & 
		(global_data['dist_to_traveler'].abs() == min_global_snps)][1].values)])
	global_output = f'{min_global_snps}{min_global_snps_label}	{min_global_snps_count}	{n_global_seqs}	{min_global_snps_locs}'
	# row output
	row_output = row[0].replace('P', '') + '\t' + row[1] + '\t' + row[2].replace('USA', ', USA') + '\t\t' + travel_output+'\t'+ga_output+'\t'+global_output+'\n'
	all_output += row_output



with open('test.tsv', 'w') as fp:
	fp.write(all_output)