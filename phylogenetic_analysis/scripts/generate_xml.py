#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd
from treetime.utils import numeric_date
import argparse


def generate_sequence_block(sequence_list):
    sequence_block = ''
    for item in sequence_list:
        sequence_block += '<sequence id="seq_{0}" spec="Sequence" taxon="{0}" \
        totalcount="4" value="{1}" />\n'.format(item.id, item.seq)
    return(sequence_block)


def generate_time_block(sequence_list, date_dict):
    time_block = ''
    for item in sequence_list:
        time_block += f'{item.description}={date_dict[item.description]},'
    # removes trailing ','
    time_block = ''.join(list(time_block)[0:-1])
    return(time_block)


def generate_trait_block(sequence_list, trait_dict):
    trait_block = ''
    for item in sequence_list:
        trait_block += f'{item.description}={trait_dict[item.description]},'
    # removes trailing ','
    trait_block = ''.join(list(trait_block)[0:-1])
    return(trait_block)


def generate_code_map(code_list):
    all_values = "? = "
    all_values += ' '.join([str(i[0]) for i in code_list])
    code_map = ''.join([str(i[1])+"="+str(i[0])+',' for i in code_list])
    # remove trailing ,
    code_map+=all_values
    return(code_map)


def generate_xml():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--xmlTemplate', 
        help='path to xml template file')
    parser.add_argument('--seqs', 
        help='path to aln')
    parser.add_argument('--seqNameSep',
        help='character to split sequence name on to match metadata',
        default='|')
    parser.add_argument('--seqNameField',
        help='field of split sequence name to match metadata',
        default=1, type=int)
    parser.add_argument('--metadata',
        help='path to metadata')
    parser.add_argument('--metadataDelim',
        help='metadata file delimiter',
        default='\t')
    parser.add_argument('--metadataNameCol',
        help='column in metadata with sequence names',
        type=int,
        default=1)
    parser.add_argument('--metadataDateCol',
        help='column in metadata with date',
        type=int,
        default=2)
    parser.add_argument('--metadataDateFmt',
        help='date format',
        default='%Y-%m-%d')
    parser.add_argument('--metadataTraitCol',
        help='column in metadata with trait info (fine scale resolution)',
        type=int,
        default=6)
    parser.add_argument('--traitName',
        help='name of discrete trait',
        default='discrete_trait')
    parser.add_argument('--alnName',
        help='name of discrete trait',
        default='sars_cov_2')
    parser.add_argument('--outName',
        help='name of output file',
        default='sars_cov_2')
    args = parser.parse_args()
    # default values used for testing
    #args.seqs = 'data/weighted_downsampling_old/19B/19B.fasta'
    #args.metadata = 'data/weighted_downsampling_old/19B/19B_seqs.tsv'
    #args.minTraitSize = 5
    #args.metadataTraitColPrimary = 7
    #args.metadataTraitColSecondary = 3
    seqs = list(SeqIO.parse(args.seqs, 'fasta'))
    n_seqs = len(seqs)
    print(f'{n_seqs} sequences in input fasta')
    metadata = \
        pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
    metadata_names = set(metadata[args.metadataNameCol])
    seqs = [i for i in seqs if i.description.split(args.seqNameSep)[args.seqNameField] in metadata_names]
    seq_names = set([i.description.split(args.seqNameSep)[args.seqNameField] for i in seqs])
    metadata[metadata[args.metadataNameCol].isin(seq_names)][
        args.metadataTraitCol].value_counts().to_csv(
            args.outName+'_trait_counts.tsv', sep='\t', header=None)
    print(f'{len(seqs)} sequences found in metadata')
    # removes "-" from seq names, these interfere for biopython nexus parsing
    for seq in seqs:
        seq.description = \
            seq.description.replace('-', '_')
        seq.id = \
            seq.id.replace('-', '_')
        seq.name = \
            seq.name.replace('-', '_')
    name_dict = \
        {i.description.split(args.seqNameSep)[args.seqNameField]: i.description for i in seqs}
    # removes "-" from metadata names, these interfere for biopython nexus parsing
    metadata[args.metadataNameCol] = \
        metadata[args.metadataNameCol].str.replace('-', '_')
    # creates a new column with trait labels 
  
    # build trait dictionary
    trait_dict = \
        {name_dict[row[args.metadataNameCol]]: row[args.metadataTraitCol] for idx, row in metadata.iterrows()}
    # build date dict
    metadata['numeric_date'] = \
        pd.to_datetime(metadata[args.metadataDateCol], format=args.metadataDateFmt).apply(numeric_date)
    date_dict = \
        {name_dict[row[args.metadataNameCol]]: row['numeric_date'] for idx, row in metadata.iterrows()}
    sequence_block = generate_sequence_block(seqs)
    time_block = generate_time_block(seqs, date_dict)
    trait_block = generate_trait_block(seqs, trait_dict)
    trait_code_map = \
        generate_code_map([(idx, i) for idx, i in enumerate(set(trait_dict.values()))])
    n_traits = int(len(set(trait_dict.values())))
    n_traits_less_1 = int(n_traits -1)
    n_trait_dimensions = int(n_traits*(n_traits-1)/2)
    with open(args.xmlTemplate, 'r') as infile:
            template = infile.read()
    template = template.replace('<!-- ALN_NAME -->', args.alnName)
    template = template.replace('<!-- TRAIT_NAME -->', args.traitName)
    template = template.replace('<!-- SEQUENCE_BLOCK -->', sequence_block)
    template = template.replace('<!-- TIME_BLOCK -->', time_block)
    template = template.replace('<!-- TRAIT_BLOCK -->', trait_block)
    template = template.replace('<!-- TRAIT_CODE_MAP -->', trait_code_map)
    template = template.replace('<!-- N_TRAITS -->', str(n_traits))
    template = template.replace('<!-- N_TRAITS_LESS_1 -->', str(n_traits_less_1))
    template = template.replace('<!-- N_TRAIT_DIMENSIONS -->', str(n_trait_dimensions))
    template = template.replace('<!-- N_TRAIT_DIMENSIONS/4 -->', str(int(n_trait_dimensions/4)))
    template = template.replace('<!-- TRAIT_FREQS -->', str(1/n_traits))
    with open(args.outName+'.xml', 'w') as outfile:
        outfile.write(template)


if __name__ == '__main__':
    generate_xml()
