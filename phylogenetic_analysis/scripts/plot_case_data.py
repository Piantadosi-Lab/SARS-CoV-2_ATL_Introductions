import matplotlib.pyplot as plt
import pandas as pd
#from scripts.utils import plot_style
from utils import plot_style
import matplotlib.dates as mdates
from datetime import datetime
import string
import argparse


def process_case_data(case_data, focal_case_state, max_date):
        dat = pd.read_csv(case_data)
        state_dat = dat[dat['Province_State']==focal_case_state].iloc[:,11:].sum(axis=0).diff()
        state_dat.index = pd.to_datetime(state_dat.index)
        state_dat = state_dat[state_dat.index <= max_date]
        state_dat = state_dat.loc[state_dat[((state_dat != 0.0) & (~pd.isnull(state_dat)))].index[0]:]
        state_dat.to_csv(f'data/time_series_{focal_case_state}.tsv', sep='\t', header=None)
        county_dat = dat[dat['Province_State']==focal_case_state].iloc[:,10:]
        county_dat = county_dat.set_index('Combined_Key')
        county_dat.columns = pd.to_datetime(county_dat.columns)
        county_dat = county_dat.loc[:,[i for i in county_dat.columns if i < datetime.strptime(max_date, '%Y-%m-%d')]]
        county_dat.sum(axis=1).sort_values().to_csv(f'data/cumulative_county_{focal_case_state}.tsv', sep='\t', header=None)
        return(state_dat, county_dat)


def process_seq_data(seq_data):
        dat = pd.read_csv(seq_data, sep='\t', header=None)
        dat[2] = pd.to_datetime(dat[2])
        seqs_per_day = dat.groupby(2).size()
        seqs_per_day.to_csv('.'.join(seq_data.split('.')[0:-1])+'_count.tsv', sep='\t', header=None)
        return(seqs_per_day)




def run():
        parser = argparse.ArgumentParser()
        parser.add_argument('--caseData')
        parser.add_argument('--focalCaseState')
        parser.add_argument('--focalSeqState')
        parser.add_argument('--focalSeqCol')
        parser.add_argument('--maxDate')
        parser.add_argument('--seqData')
        args = parser.parse_args()
        #args.caseData = 'data/time_series_covid19_confirmed_US.csv'
        #args.seqData = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv'
        args.focalCaseState = 'Georgia'
        args.maxDate = '2020-03-31'
        plot_style()
        state_dat, county_dat = process_case_data(args.caseData, args.focalCaseState, args.maxDate)
        seqs_per_day = process_seq_data(args.seqData)

        early_date = min(seqs_per_day.index.min(), state_dat.index.min())
        late_date = max(seqs_per_day.index.max(), state_dat.index.max())
        idx = pd.date_range(early_date, late_date)

        seqs_per_day = seqs_per_day.reindex(idx, fill_value=0)
        state_dat = state_dat.reindex(idx, fill_value=0)

        fig, ax = plt.subplots(figsize=(6.4*1.25, 4.8), constrained_layout=True)
        ax2 = ax.twinx()
        ax.plot(state_dat.index, state_dat, color='#8B3C37')
        ax2.plot(seqs_per_day.index, seqs_per_day, color='#37868B')
        max_y_ticks = int(ax2.get_yticks()[-1])
        #ax2.set_yticks(range(0, max_y_ticks, 2))
        ax.set_ylabel('Cases/day', color='#8B3C37')
        ax2.set_ylabel('Sequences/day', color='#37868B')
        ax.set_xlabel('Date (2020)')
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
        ax.axvline(datetime.strptime('2020-03-15', '%Y-%m-%d'), ls='--', color='#d8dee9')
        ax.set_ylim(0, 1200)
        ax2.set_ylim(0, 120)
        for ax_idx, ax in enumerate([ax]):
                ax.text(-0.15, 1.05, string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
                        size=20, weight='bold', va="top")

        fig.savefig('figures/figure_1a.pdf')
        plt.close()



if __name__ == "__main__":
        run()


