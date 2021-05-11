import argparse
import pandas as pd
import numpy as np
import glob


def process_results(file, results):
    file_name = file.split('/')[-1].split('.')[0]
    dat = pd.read_csv(file, index_col = 0)
    dat = pd.melt(dat, id_vars=['destination region', 'source region', 'descendant tips', 'source node', '# descendants'], 
            value_vars=['midroot time', 'lower midroot time', 'upper midroot time'])
    dat['file'] = file_name
    dat = dat.sort_values('value')
    dat = dat.assign(cumsum=dat.groupby(['destination region', 'variable']).cumcount()+1)
    results = results.append(dat)
    return(results.reset_index(drop=True))


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
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--treeBootstrapDir', 
                        help='directory with bootstrap replicate trees',
                        default=None)
    args = parser.parse_args()
    #args.treeBootstrapDir = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv'
    results = pd.DataFrame()
    for file in glob.glob(args.treeBootstrapDir):
        results = process_results(file, results)
    
    midroot_introductions = \
        results[results['variable'] == 'midroot time'].groupby(
            ['file', 'destination region', 'variable'])['cumsum'].agg('max').reset_index()
    n_introductions = hpd(midroot_introductions['cumsum'], [0.025, 0.5, 0.975])
    print(n_introductions)






if __name__ == "__main__":
    run()








