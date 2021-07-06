import argparse
import numpy as np
import pandas as pd
import csv

# mostly from arviz
# but works directly with a pandas object
# and also returns median
# https://github.com/arviz-devs/arviz/blob/master/arviz/stats/stats.py
def hpd(dat, qs):
    width = qs[2] - qs[0]
    # sorts from smallest to largest
    dat = dat.sort_values()
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
    interval_width = dat[interval_idx_inc:].reset_index(drop=True) - \
                     dat[:n_intervals].reset_index(drop=True)
    # find the shortest interval
    min_idx = interval_width.values.argmin()
    hpd_interval = (dat.iloc[min_idx], dat.iloc[min_idx+interval_idx_inc])
    dat_hpd = dat[dat.between(hpd_interval[0],
                              hpd_interval[1])]
    dat_mid = dat_hpd.quantile(qs[1])
    return((hpd_interval[0], dat_mid, hpd_interval[1]))


def make_table(rows, data, qs, out_path):
    w = str(int((qs[2] - qs[0])*100))
    table = [[rows[row]['label'], rows[row]['prior']] for row in rows.keys()]
    table.insert(0, ['Parameter', 'Prior', f'Estimate ({w}% HPD)']) 
    for i, key in enumerate(rows.keys()):
        dat = data[f'{key}_HPD']
        if rows[key]['sci']:
            median = "{:.{}e}".format(dat.median(), rows[key]['decimals'])
            low_limit = "{:.{}e}".format(min(dat), rows[key]['decimals'])
            high_limit = "{:.{}e}".format(max(dat), rows[key]['decimals'])
        else:
            median = "{:.{}f}".format(dat.median(), rows[key]['decimals'])
            low_limit = "{:.{}f}".format(min(dat), rows[key]['decimals'])
            high_limit = "{:.{}f}".format(max(dat), rows[key]['decimals'])
        table[i+1].append(f'{median} ({low_limit} {high_limit})')
        with open(out_path, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerows(table)


def process_results(file_path, burnin, qs):
    results = {}
    log = pd.read_csv(file_path, sep="\t", comment="#")
    #plot_trace(log, args)
    log = log[log['Sample']>=max(log['Sample'])*burnin]
    for d in log.columns[1:]:   # excludes sample column
        results[d] = log[d]
        limits = hpd(log[d], qs)
        results[d+'_HPD'] = log[d][log[d].between(limits[0], 
                                                     limits[2])]
    return(results)


             

def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--log', 
        help='beast log file to parse',
        default=None)
    parser.add_argument('--burnin', 
        help='what % to discard as burnin',
        default=0.10, type=float)
    parser.add_argument('--qs', 
        help='what quantiles to report',
        default=[0.025, 0.5, 0.975], type=float, nargs=3)
    args = parser.parse_args()
    #args.log = 'data/19B_subclade/19B_subclade.log'
    table_rows = {  'posterior':    
                        {'label':           'Posterior',
                         'decimals':        0, 
                         'sci':             False,
                         'prior':           ''}, 
                    'TreeHeight':   
                        {'label':           'Tree height',
                         'decimals':        3,
                         'sci':             False,
                         'prior':           ''},
                    'proportionInvariant.19B_sublcade':   
                        {'label':           'Proportion invariant',
                         'decimals':        3,
                         'sci':             False,
                         'prior':           'Uniform(0, 1)'},
                    'ucldMean.19B_sublcade':    
                        {'label':           'UCLD Mean',
                         'decimals':        2,
                         'sci':             True,
                         'prior':           f'Normal(1.0e-3, 1.0e-4)'},
                    'ucldStdev.19B_sublcade':    
                        {'label':           'UCLD $\N{GREEK SMALL LETTER SIGMA}',
                         'decimals':        2,
                         'sci':             False,
                         'prior':           f'Exponential(0.33)'},
                    'kappa.19B_sublcade':        
                        {'label':           f'\N{GREEK SMALL LETTER KAPPA}',
                         'decimals':        3,
                         'sci':             False,
                         'prior':           'Lognormal(1.0, 1.25)'},
                    'gammaShape.19B_sublcade': 
                        {'label':           f'\N{GREEK SMALL LETTER GAMMA}',
                         'decimals':        2,
                         'sci':             True,
                         'prior':           'Exponential(1.0)'}}
                    
    all_results = process_results(args.log, args.burnin, args.qs)
    make_table(table_rows, all_results, args.qs, f'tables/{args.log.split("/")[-1].split(".")[0]}.csv')





if __name__ == "__main__":
    run()

    
