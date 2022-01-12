import argparse
import pandas as pd
from collections import Counter
import numpy as np
import datetime



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--metadata')
	parser.add_argument('--delim', default='\t')
	parser.add_argument('--dateCol', default=1, type=int)
	parser.add_argument('--groupCol', default=0, type=int)
	parser.add_argument('--timeGroup', default='week')
	parser.add_argument('--nPerGroup', type=int, default=1)
	parser.add_argument('--nReps', type=int, default=1)
	args = parser.parse_args()
	#args.metadata = 'data/19B_subclade/19B_subclade_beast_include.tsv'
	#args.nPerGroup = 2
	#args.groupCol = 7
	#args.dateCol = 2
	dat = pd.read_csv(args.metadata, sep=args.delim, header=None)
	print(dat)
	dat[args.dateCol] = pd.to_datetime(dat[args.dateCol])
	if args.timeGroup.lower() == 'day':
		dat['end_date'] = pd.to_datetime(dat[args.dateCol])
	elif args.timeGroup.lower() == 'week':
		dat['end_date'] = \
			  pd.to_datetime(dat[args.dateCol]).apply(lambda k:
			  k+datetime.timedelta(days= 6 - k.weekday()))
	elif args.timeGroup.lower() == 'month':
		dat['end_date'] = \
			  pd.to_datetime(dat[args.dateCol]).apply(lambda k: str(k.year) + '-' + str(k.month))
	elif args.timeGroup.lower() == 'year':
		dat['end_date'] = \
			  pd.to_datetime(dat[args.dateCol]).apply(lambda k: k.year)
	else:
		dat['end_date'] = 0
	n_take = {group: min(group_dat.shape[0], args.nPerGroup) for 
		group, group_dat in dat.groupby([args.groupCol, 'end_date'])}
	for rep in range(args.nReps):
		np.random.seed(rep)
		sampled_dat = dat.groupby([args.groupCol, 'end_date']).apply(lambda x: 
			x.iloc[np.random.choice(x.shape[0], 
				size=n_take[x.name], replace=False)]).reset_index(drop=True)
		sampled_dat.to_csv('.'.join(
			args.metadata.split('.')[:-1])+f'_sampled_{args.nPerGroup}_rep_{rep}.tsv',
			sep='\t', header=None, index=None)
        






if __name__ == "__main__":
    run()
