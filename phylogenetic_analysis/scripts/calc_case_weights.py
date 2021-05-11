import pandas as pd
import datetime
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
import os


def process_dat(dat, groupby_var, max_date, datetime_fmt):
    # groups and sums
    dat = dat.groupby(groupby_var).agg('sum')
    # converts columns to datetime
    dat.columns =  pd.to_datetime(dat.columns, format=datetime_fmt)
    # gets only through march 31st
    dat = \
        dat.loc[:,dat.columns <= 
            datetime.datetime.strptime(max_date, datetime_fmt)].\
        sort_values(max_date, ascending=False)
    top_dat = dat.loc[dat.index[0:9],:]
    summarized_dat = pd.concat([top_dat, 
        pd.DataFrame(dat.sum() - top_dat.sum(), 
            columns=['Other']).transpose()])
    return(dat, summarized_dat)


def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--globalDat', 
                        help='csv file with country cases by day (from hopkins CCSE)',
                        default=None)
    parser.add_argument('--usDat', 
                        help='csv file with US cases by day (from hopkins CCSE)',
                        default=None)
    parser.add_argument('--dateMax', 
                        help='max date',
                        default=None)
    parser.add_argument('--dateFmt', 
                        help='date format',
                        default='%m/%d/%y')
    parser.add_argument('--focalStates',
        help='focal states, exlcuded from USA weight',
        nargs='+')
    args = parser.parse_args()
    # /// PROCESS GLOBAL DAT ///
    #args.dateMax='03/31/20'
    #args.globalDat = 'data/time_series_covid19_confirmed_global.csv'
    #args.usDat = 'data/time_series_covid19_confirmed_US.csv'
    #args.focalStates = ['Georgia']
    global_dat = pd.read_csv(args.globalDat, sep=',')
    # to make sure hopkins definitions line up with gisaid
    # some hopkins provinces are broken out as GISAID countries
    province_country_mappings = {'Curacao': 'Curacao', 
        'Faroe Islands': 'Faroe Islands',
        'Aruba': 'Aruba', 
        'St Martin': 'Saint Martin', 
        'Hong Kong': 'Hong Kong', 
        'Saint Barthelemy': 'Saint Barthélemy',
        'Guadeloupe': 'Guadeloupe',
        'Reunion': 'Reunion',
        'Bermuda': 'Bermuda',
        'St Martin': 'St Martin',
        'Sint Maarten': 'Sint Maarten'}
    global_dat.loc[global_dat['Province/State'].isin(province_country_mappings.keys()), 
            'Country/Region'] = \
        global_dat.loc[global_dat['Province/State'].isin(province_country_mappings.keys()), 
            'Province/State'].map(province_country_mappings)
    # some countries are spelled different in GISAID
    country_mappings = {'Korea, South': 'South Korea', 
        'Congo (Brazzaville)': 'Republic of the Congo',
        'Congo (Kinshasa)': 'Democratic Republic of the Congo',
        'Czechia': 'Czech Republic', 
        }
    global_dat.loc[global_dat['Country/Region'].isin(country_mappings), 
            'Country/Region'] = \
        global_dat.loc[global_dat['Country/Region'].isin(country_mappings), 
            'Country/Region'].map(country_mappings)
    # removes unneeded columns
    global_dat = global_dat[global_dat.columns[0:2].tolist() + 
        global_dat.columns[4:].tolist()]
    all_country_dat, summarized_country_dat = \
        process_dat(global_dat, 'Country/Region', args.dateMax, args.dateFmt)
    # finally, add to dictionary
    country_case_dict = {idx: i[-1] for idx, i in all_country_dat.iterrows()}
    # some countries are spelled different in GISAID
    # and hopkins
    # will not double count because these are all
    # strict 1-1 mappings
    # (confirmed by hand)
    collapse_countries = {
        "Côte d'Ivoire": "Cote d'Ivoire",
        'Myanmar': 'Burma',
        'Saint Martin': 'St Martin'}
    for key, value in collapse_countries.items():
        country_case_dict[key] = all_country_dat.loc[value,:][-1]
    # //// PROCESS US STATE DAT ////
    # gets US data
    us_dat = pd.read_csv(args.usDat, sep=',')
    # removes unneeded columns
    us_dat = us_dat[['Province_State'] + us_dat.columns[11:].tolist()]
    all_state_dat, summarized_state_dat = \
        process_dat(us_dat, 'Province_State', args.dateMax, args.dateFmt)
    # //// PROCESS SEQUENCE DAT ////
    # some US territories are listed as part of 
    # the US in the hopkins data but broken out as country
    # in GISAID
    # so we want to dictionary and remove from us state data
    for country in ['Guam']:
        # already filtered by max date so can take -1 index
        country_dat = all_state_dat.loc[country,:].iloc[-1]
        country_case_dict[country] = country_dat
        all_state_dat = all_state_dat[all_state_dat.index != country]
    # remove focal states -- likely these are force included
    # and thus do not want to count their cases in the weights
    for state in args.focalStates:
        all_state_dat = all_state_dat[all_state_dat.index != state]
    # then, sum over remaining US states to get USA weights
    country_case_dict['USA'] = all_state_dat.sum()[-1]
    # finally, we convert to DF, calc relative weights, and save
    country_case_dat = pd.DataFrame(zip(country_case_dict.keys(), country_case_dict.values()))
    country_case_dat[2] = country_case_dat[1]/country_case_dat[1].sum()
    country_case_dat.to_csv('data/country_case_weights.tsv', header=None, index=None, sep='\t')



if __name__ == "__main__":
    run()

    
'''
all_state_dat_summed = all_state_dat.sum()

    
    # finally, calcualte sequence weights
    country_case_dict = {idx: i[-1] for idx, i in all_country_dat.iterrows()}
    state_case_dict = {'USA_'+idx: i[-1] for idx, i in all_state_dat.iterrows()}
    region_case_dict = {**country_case_dict, **state_case_dict}
    sequence_cases = \
        pd.DataFrame.from_dict({key: region_case_dict[key] for key in sequence_locs if key in region_case_dict.keys()}, 
            orient='index').sort_values(by=0, ascending=False)
    sequence_cases['weight'] = \
        sequence_cases[0]/sum(sequence_cases[0])
    #sequence_cases['weight_log'] = \
    #    np.log10(sequence_cases[0]+1)/sum(np.log10(sequence_cases[0]+1))
    sequence_cases.to_csv(os.path.splitext(args.sequences)[0]+'_weights.csv', header=None)
    # get
    # plots country dat
    fig, axs = plt.subplots(1, 2, figsize=(6.4*2,4.8), constrained_layout=True)
    summarized_country_dat.transpose().plot.area(ax=axs[0], 
        colormap=sns.color_palette("tab20", as_cmap=True), alpha=0.5, lw=1)
    axs[0].set_ylabel('cumulative cases')
    axs[0].set_title('global')
    summarized_state_dat.transpose().plot.area(ax=axs[1], 
        colormap=sns.color_palette("tab20b", as_cmap=True), alpha=0.5, lw=1)
    axs[1].set_title('usa')
    fig.savefig('cumulative_cases.pdf')
    plt.close()

# //// NOT OVERSAMPLING US SEQUENCES ////
args.globalDat = 'data/time_series_covid19_confirmed_global.csv'
args.metadata = 'data/gisaid_hcov-19_2021_01_26_prop_date_location_coarse.csv'
args.dateFmt = '%m/%d/%y'
args.dateMax = '03/31/20'
global_dat = pd.read_csv(args.globalDat, sep=',')
# to make sure hopkins definitions line up with gisaid
# some provinces are broken out as GISAID countries
province_country_mappings = {'Curacao': 'Curacao', 
    'Faroe Islands': 'Faroe Islands',
    'Aruba': 'Aruba', 
    'St Martin': 'Saint Martin', 
    'Hong Kong': 'Hong Kong', 
    'Saint Barthelemy': 'Saint Barthélemy',
    'Guadeloupe': 'Guadeloupe',
    'Reunion': 'Reunion'}
global_dat.loc[global_dat['Province/State'].isin(province_country_mappings.keys()), 
        'Country/Region'] = \
    global_dat.loc[global_dat['Province/State'].isin(province_country_mappings.keys()), 
        'Province/State'].map(province_country_mappings)
# some countries are spelled different in GISAID
country_mappings = {'Korea, South': 'South Korea', 
    'Congo (Brazzaville)': 'Republic of Congo',
    'Congo (Kinshasa)': 'Democratic Republic of the Congo',
    'Czechia': 'Czech Republic'}
global_dat.loc[global_dat['Country/Region'].isin(country_mappings), 
        'Country/Region'] = \
    global_dat.loc[global_dat['Country/Region'].isin(country_mappings), 
        'Country/Region'].map(country_mappings)
# removes unneeded columns
global_dat = global_dat[global_dat.columns[0:2].tolist() + 
    global_dat.columns[4:].tolist()]
all_country_dat, summarized_country_dat = \
    process_dat(global_dat, 'Country/Region', args.dateMax, args.dateFmt)


# //// PROCESS US STATE DAT ////
args.usDat = 'data/time_series_covid19_confirmed_US.csv'
# gets US data
us_dat = pd.read_csv(args.usDat, sep=',')
# removes unneeded columns
us_dat = us_dat[['Province_State'] + us_dat.columns[11:].tolist()]
all_state_dat, summarized_state_dat = \
    process_dat(us_dat, 'Province_State', args.dateMax, args.dateFmt)
non_ga_dat = all_state_dat.sum() - all_state_dat.loc['Georgia']

country_case_dict = {idx: i[-1] for idx, i in all_country_dat.iterrows()}
# adds US data
country_case_dict.update({'Other(USA)': non_ga_dat[-1], 
    'Georgia(USA)': all_state_dat.loc['Georgia'][-1]})

metadata = pd.read_csv(args.metadata, sep=',', header=None)
collapse_countries = {'Palestine': 'Israel', 
    'Viet Nam': 'Vietnam',
    'Czech Republic': 'Czech republic',
    'Crimea': 'Russia'}
sequence_locs = \
    [collapse_countries[i] if i in collapse_countries.keys() else i for i in metadata[2].unique()]
sequence_cases = \
    pd.DataFrame.from_dict({key: country_case_dict[key] for key in sequence_locs if key in country_case_dict.keys()}, 
        orient='index').sort_values(by=0, ascending=False)

sequence_cases['weight'] = \
    sequence_cases[0]/sum(sequence_cases[0])
sequence_cases.to_csv(os.path.splitext(args.metadata)[0]+'_coarse_weights.csv', header=None)



# plots sampled sequences against number of cases
sequences = list(SeqIO.parse('data/ga_focused_prop_aligned_ref_filtered_masked.fasta', 'fasta'))
metadata = pd.read_csv('data/ga_focused_prop_date_location_name.csv', sep=',', header=None)
case_counts = pd.read_csv('data/gisaid_hcov-19_2021_01_26_weights.csv', sep=',', header=None)
seq_names = pd.DataFrame([i.description.split('|')[1] for i in sequences])
seq_counts = \
    seq_names.merge(metadata, left_on=0, right_on=0).groupby(2).apply(lambda x: len(x)).reset_index()

seq_counts['prop_seqs'] = seq_counts[0]/seq_counts[0].sum()
seq_data = seq_counts.merge(case_counts, left_on=2, right_on=0)

fig, axs = plt.subplots(1,1, figsize=(6.4, 4.8), constrained_layout=True)
axs.scatter(seq_data[1], seq_data['0_x'], color='#333333', alpha=0.85)
axs.set_xlabel('cumulative cases')
axs.set_ylabel('number of down sampled sequences')
fig.savefig('figures/sampled_sequences.pdf')

seq_dates = \
    seq_names.merge(metadata, left_on=0, right_on=0).groupby(1).apply(lambda x: len(x)).reset_index()

seq_dates[1] = pd.to_datetime(seq_dates[1])


fig, axs = plt.subplots(constrained_layout=True)
axs.scatter(seq_dates[1], seq_dates[0],
    color='#333333', alpha=0.85)
axs.set_xlabel('date')
axs.set_ylabel('numer of down sampled sequences')
axs.tick_params(labelrotation=45, axis='x')
fig.savefig('figures/sampled_sequence_dates.pdf')
'''

