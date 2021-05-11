import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import datetime
from utils import convert_to_datetime
#from scripts.utils import convert_to_datetime
from treetime.utils import numeric_date
import sys 
from collections import Counter


def run():
    parser = argparse.ArgumentParser()
    # input files
    # //// FASTA FILE PARSING ////
    parser.add_argument('--sequences', 
        default=None)
    parser.add_argument('--seqNameSep',
        default='|')
    parser.add_argument('--seqNameField',
        default=1, type=int, 
        help='which field in the sequence name corresponds to the metadata name column')
    # //// METADATA FILE PARSING ////
    parser.add_argument('--metadata', 
        default=None, help='metadata file path, assumes no header')
    parser.add_argument('--metadataDelim',
        default='\t', 
        help='metadata delimiter')
    parser.add_argument('--metadataIDCol',
        default=1,
        type=int,
        help='which column in the metadata file contains sequence names (corresponds to seqNameField)')
    parser.add_argument('--metadataDateCol',
        default=2,
        type=int,
        help='which column in the metadata file contains sequece dates')
    parser.add_argument('--metadataDateFmt',
        default='%Y-%m-%d',
        help='format of the date column')
    parser.add_argument('--metadataRegionCol',
        help='which column in the metadata file contains the region (must correspond with region weights)',
        default=4,
        type=int)
    # //// DOWNSAMPLE PARAMETERS ////
    parser.add_argument('--minDate', 
        default=None)
    parser.add_argument('--maxDate', 
        default=None)
    parser.add_argument('--samplesPerWeek', 
        default=None, type=int)
    parser.add_argument('--regionWeights', 
        default=None)
    parser.add_argument('--targetN', 
        default=100, type=int, help='total number of sequences in alignment to target when using weights')
    parser.add_argument('--interRegionWeights', 
        default=None, help='weights to use to sample sequences within each region')
    parser.add_argument('--exclude', 
        help='which countries to exclude, here overrides include',
        default=None)
    parser.add_argument('--include', 
        help='which items to include in the downsampling',
        default=None)
    parser.add_argument('--goodChars', 
        help='which characters to allow in fasta file',
        default=['A', 'C', 'G', 'T', '-', 'N', 'R', 'Y', 'S', 'W', 
        'K', 'M', 'D', 'H', 'B', 'V', '?'],
        nargs='+')
    parser.add_argument('--outName', 
        help='outname to save data to if None, defaults to modified input fasta name',
        default=None)
    args = parser.parse_args()
    #args.sequences = 'data/gisaid_hcov-19_2021_01_26_aligned_ref_filtered_masked.fasta'
    #args.metadata = './data/metadata_aligned.tsv'
    #args.include = './config/include.tsv'
    #args.exclude = './config/exclude.csv'
    #args.maxDate = '03-31-2020'
    #args.regionWeights = 'data/country_case_weights.tsv'
    #args.targetN = 6000
    #args.interRegionWeights = 'data/gisaid_hcov-19_2021_01_26_EHC_aligned_ref_filtered_masked_min_dist_weights.tsv'
        # read in metadata
    metadata=pd.read_csv(args.metadata, 
        sep=args.metadataDelim, header=None, low_memory=False)
    # read in sequences 
    seqs = list(SeqIO.parse(args.sequences, 'fasta'))
    # //// 1) removes sequences from sequence file with bad characters ////
    # todo faster
    # todo move up
    good_chars_str = ''.join(set(args.goodChars))
    len_seqs = len(seqs)
    seqs = [i for i in seqs if not i.seq.strip(good_chars_str)]
    print(f'{len_seqs - len(seqs)} sequences removed for non nucleotide characters', 
        file=sys.stderr)
    seqs_names = set([i.description.split(args.seqNameSep)[args.seqNameField] for i in seqs])
    
    # //// 2) remove metadata entires not in the fasta file ////
    len_metadata = metadata.shape[0]
    metadata = metadata[metadata[args.metadataIDCol].isin(seqs_names)]
    print(f'{len_metadata - metadata.shape[0]} metadata entries removed because they are not in the fasta file', 
        file=sys.stderr)    

    # //// 3) remove metadata entries in the manually excluded categories ////
    len_metadata = metadata.shape[0]
    if args.exclude:
        to_exclude = pd.read_csv(args.exclude, sep=args.metadataDelim, header=None, comment='#')
        for idx,row in to_exclude.iterrows():
            metadata = \
                metadata[metadata[row[0]].str.contains(row[1]) == False]
        print(f'{len_metadata - metadata.shape[0]} sequences excluded by exclude file', 
            file=sys.stderr)
        len_metadata = metadata.shape[0]

    # //// 4) convert metadata date columns to datetime, remove ambiguous dates, removes anything outside date range ////
    # converts date column to datetime
    metadata[args.metadataDateCol] = \
        metadata[args.metadataDateCol].apply(lambda k: convert_to_datetime(k, args.metadataDateFmt))
    # filters ambiguous dates
    metadata = metadata[metadata[args.metadataDateCol].notnull()]
    print(f'{len_metadata - metadata.shape[0]} sequences removed due to ambiguous dates', 
        file=sys.stderr)
    len_metadata = metadata.shape[0]
    # removes anything dated after max_date or before min_date assuming they're given
    if args.maxDate:
        metadata = \
            metadata[metadata[args.metadataDateCol] <= args.maxDate]
    if args.minDate: 
        metadata = \
            metadata[metadata[args.metadataDateCol] >= args.minDate]
    print(f'{len_metadata - metadata.shape[0]} sequences removed by date restriction', file=sys.stderr)

    # //// 5) get sequences in manually included list, only if they have non-ambiguous dates ////
    # get manually included sequences
    sampled_metadata = pd.DataFrame()
    if args.include:
        to_include = \
            pd.read_csv(args.include, sep=args.metadataDelim, header=None, comment='#')
        for idx,row in to_include.iterrows():
            row_metadata = metadata[metadata[row[0]] == row[1]]
            sampled_metadata = pd.concat([sampled_metadata, 
                row_metadata])
            # removes the included rows from the metadata to avoid double counting
            metadata = metadata.drop(row_metadata.index)
        #for idx,col in to_include.groupby(0):
        #    include_items = set(col[1])            
        #    sampled_metadata = pd.concat([sampled_metadata,
        #        metadata[metadata[].isin(include_items)]])
        #    metadata = metadata[metadata[idx].isin(include_items) == False]
        print(f'{sampled_metadata.shape[0]} sequences added by include file', file=sys.stderr)
    
    # //// 6) finally, downsample metadata ///
    # args.SamplesPerWeek and args.weights are currently incompatible
    if args.samplesPerWeek and args.regionWeights:
        raise exception('samplesPerWeek and weights are currently incompatible')
    elif args.regionWeights:
        region_weights = pd.read_csv(args.regionWeights, sep=args.metadataDelim, header=None)
        # convert weights to target # sequences per region given target size
        region_weights['target_n_seqs'] = region_weights[2]*args.targetN
        # get number of sequences belonging to each location in sequence file
        n_avail_sequences = Counter(metadata[args.metadataRegionCol])
        region_weights['n_avail'] = region_weights[0].map(n_avail_sequences)
        region_weights['take'] = \
            region_weights[['target_n_seqs', 'n_avail']].min(axis=1).astype(int)
        n_choose_dict = {i[0]: i['take'] for idx, i in region_weights.iterrows()}
        # in addition to weighting regions, we can weight the sequences within each region
        if args.interRegionWeights:
            inter_region_weights = pd.read_csv(args.interRegionWeights, sep=args.metadataDelim, header=None)
            inter_region_weights_dict = {i[0]: i[1] for idx, i in inter_region_weights.iterrows()}
            metadata['weight'] = metadata[args.metadataIDCol].map(inter_region_weights_dict)
        else:
            metadata['weight'] = 1
        # Palestine broken out from Israel in GISAID but not case data
        # thus, we need to combine for downsampling purposes only
        metadata['region_downsampling'] = metadata[args.metadataRegionCol]
        metadata.loc[metadata['region_downsampling'] == \
            'Palestine','region_downsampling'] = 'Israel'
        weight_sampled_metadata =\
             metadata.groupby('region_downsampling').apply(lambda x:
                x.iloc[np.random.choice(np.arange(len(x)),
            size=n_choose_dict[x.name] if x.name in n_choose_dict.keys() else 0, 
            replace=False, p=x['weight']/x['weight'].sum())])
        sampled_metadata = pd.concat([sampled_metadata,
           weight_sampled_metadata])
    elif args.samplesPerWeek:
        # if args.SamplesPerWeek assume downsampling is based on n samples per country per week ////
        # takes a set of names so as not to double count the manually included regions
        # todo make this better, remove manualyl included entries from metadata
        metadata['week'] = metadata[args.metadataDateCol].dt.isocalendar().week
        sampled_metadata = pd.concat([sampled_metadata,\
            metadata.groupby([args.metadataRegionCol, 'week']).apply(
            lambda x: x.iloc[np.random.choice(np.arange(len(x)), 
                size=min(len(x), args.samplesPerWeek), replace=False)])])
    sampled_seq_names = set(sampled_metadata[args.metadataIDCol])
    print(f'{len(sampled_seq_names)} samples selected for downsampling', file=sys.stderr)
    # want to remove any sequences that are duplicated in the fasta
    sampled_seqs = []
    included_seqs_descriptions = set([])
    for idx, i in enumerate(seqs):
        if i.description.split(args.seqNameSep)[args.seqNameField] in sampled_seq_names:
            if i.description.split(args.seqNameSep)[args.seqNameField] not in included_seqs_descriptions:
                sampled_seqs.append(i)
                included_seqs_descriptions.add(i.description.split(args.seqNameSep)[args.seqNameField])
    print(f'{len(sampled_seq_names) - len(sampled_seqs)} samples missing in fasta file', file=sys.stderr)
    print(f'{len(sampled_seqs)} samples selected by downsampling', file=sys.stderr)
    # remove spaces from sequence names for IQ-Tree
    # todo confirm this works
    for i in sampled_seqs:
        i.description = i.description.replace(' ', '_')
        i.name = i.description
        i.id = i.description
    if args.outName: 
        print(args.outName+'.fasta')
        with open(args.outName+'.fasta', 'w') as out:
            SeqIO.write(sampled_seqs, out, 'fasta')
    else:
        print(args.sequences.replace('.fasta', f'_{len(sampled_seqs)}.fasta'), file=sys.stdout)
        with open(args.sequences.replace('.fasta', f'_{len(sampled_seqs)}.fasta'), 'w') as out:
            SeqIO.write(sampled_seqs, out, 'fasta')



if __name__ == "__main__":
    run()
