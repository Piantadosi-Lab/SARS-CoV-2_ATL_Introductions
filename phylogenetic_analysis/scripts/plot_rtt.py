import argparse
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib as mpl
from treetime.utils import numeric_date
from scipy.stats import linregress
from datetime import datetime
from utils import plot_style

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--data', 
	                    default=None)
	parser.add_argument('--outName', 
	                    default='rtt')
	args = parser.parse_args()
	#args.data = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_rtt.csv'
	plot_style()

	dat = pd.read_csv(args.data, sep=',', header=None)
	print(dat)
	regress = linregress(dat[0], dat[1])


	x_ticks = [datetime.strptime('2020-01-01', '%Y-%m-%d'),
		datetime.strptime('2020-01-15', '%Y-%m-%d'),
		datetime.strptime('2020-02-01', '%Y-%m-%d'),
		datetime.strptime('2020-02-15', '%Y-%m-%d'),
		datetime.strptime('2020-03-01', '%Y-%m-%d'),
		datetime.strptime('2020-03-15', '%Y-%m-%d'),
		datetime.strptime('2020-04-01', '%Y-%m-%d')]
	x_labels = \
		[i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
	x_ticks = [numeric_date(i) for i in x_ticks]

	fig, ax = plt.subplots(constrained_layout=True)
	ax.scatter(dat[0], dat[1], color='#333333', alpha=0.75)
	ax.plot(dat[0], regress.slope*dat[0]+regress.intercept, 
		color='#BF616A', lw=2)
	ax.set_xlabel('Date (2020)')
	ax.set_ylabel('Divergence (subs/site)')
	ax.ticklabel_format(useOffset=False)
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_labels)
	ax.set_xlim(2019.9660612695561, 2020.25+0.034)
	[ax.spines[loc].set_visible(False) for loc in ['top', 'right']]
	if regress.pvalue < 0.0001:
		ax.text(0.05, 0.85, "y = {:.2e}x - {:.2f}\n$p < 0.0001$".format(regress.slope, -1*regress.intercept), 
		transform=ax.transAxes, color='#BF616A', size=14)
	else:
		ax.text(0.05, 0.85, "y = {:.2e}*x - {:.2f}\n$p$ = {:.2e}".format(regress.slope, -1*regress.intercept, regress.pvalue), 
			transform=ax.transAxes, color='#BF616A', size=14)
	fig.savefig(f'{args.outName}.pdf')





if __name__ == "__main__":
    run()





