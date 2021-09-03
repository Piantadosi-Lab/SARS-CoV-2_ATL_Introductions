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
from utils import plot_style


def plot_table(ax, table_file=None, plot_config=None):
	x_vals = [0, 0.4, 0.7]
	y_vals = [0.75, 0.69]
	table_df = pd.read_csv(table_file, sep='\t', header=None)
	font_size = 10+int(ax.figure.bbox.height*0.01)*0.25
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
	match_seqs['week_end_date'] = \
		pd.to_datetime(match_seqs[date_col]).apply(lambda k: 
			k+datetime.timedelta(days= 6 - k.weekday()))
	n_per_week = \
		match_seqs.groupby([loc_col, 'week_end_date']).size().reset_index()
	# sums by geo and gets top N geos
	top_n_locs = \
		set(n_per_week[[loc_col, 0]].groupby([loc_col]).\
			sum().reset_index().sort_values(by=0,ascending=False)[:n_locs][loc_col])
	# add new location column
	new_location = \
		match_seqs.apply(lambda k: k[loc_col] if k[loc_col] in 
			top_n_locs else 'Other (USA)' if 'USA' in k[loc_col] else 'Other (intl.)', 
			axis=1).copy(deep=True)
	match_seqs['new_location'] = new_location
	n_per_week = \
		match_seqs.groupby(['week_end_date', 'new_location']).size().reset_index()
	n_per_week = n_per_week.pivot(index='week_end_date', columns='new_location', values=0)
	n_per_week.to_csv('.'.join(match_profile_file.split('.')[0:-1])+'_n_loc_per_week.csv')
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
	font_size = 10+int(ax.figure.bbox.height*0.01)*3.5
	ax_height = (ax.get_position().y1 - ax.get_position().y0) * ax.figure.bbox.height
	for col_idx, col in enumerate(col_order[::-1]):
		ax.text(0.69, 0.9-(1.35*font_size/ax_height)*col_idx, 
			col, color=plot_config['colors'][col], 
			transform=ax.transAxes,
			size=font_size, 
			path_effects=plot_config['effects'])
	'''
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
	'''
	ax.get_legend().remove()
	x_ticks = [
		datetime.datetime.strptime('2020-03-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-03-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-04-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-04-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-05-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-05-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-06-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-06-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-07-01', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-07-15', '%Y-%m-%d'),
		datetime.datetime.strptime('2020-08-01', '%Y-%m-%d')]
	x_labels = \
        [i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_labels, rotation=0, ha='center')
	ax.grid(axis='x', ls='-', color='#d8dee9', zorder=0)
	ax.set_axisbelow(True)
	'''
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
	'''
	return(ax)



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--mutationalProfile', 
						help='file with the mutational profile of sequences in tree')
	parser.add_argument('--matchMutationalProfile',
		help='estimated number of cases')
	parser.add_argument('--matchLocCol',
		help='which column has country data', 
		type=int, 
		default=8)
	parser.add_argument('--matchDateCol',
		help='which column has week numbers',
		type=int,
		default=3)
	parser.add_argument('--config',
		help='json file with config information (labels, colors, etc.)')
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


	fig = plt.figure(figsize=(6.4*1.5,4.8))
	#plt.rcParams['hatch.linewidth'] = 0.5
	gs=GridSpec(1,10, figure=fig)
	ax0 = fig.add_subplot(gs[:2])
	ax1 = fig.add_subplot(gs[2:])

	ax0 = plot_table(ax0, table_file=args.mutationalProfile, plot_config=plot_config)
	ax1 = \
		plot_count_seqs_time(ax1, match_profile_file=args.matchMutationalProfile, 
			loc_col=args.matchLocCol, date_col=args.matchDateCol, 
			n_locs=6,  plot_config=plot_config)
	y_poses = [1.0, 1.0, 1.1]
	for ax_idx, ax in enumerate([ax0, ax1]):
	    ax.text(-0.09, y_poses[ax_idx], string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
	            size=18, weight='bold', va="top")
	    ax.set_axisbelow(True)
	fig.tight_layout()
	fig.savefig(f'{plot_config["out_name"]}.pdf')
	plt.close()


if __name__ == "__main__":
    run()
