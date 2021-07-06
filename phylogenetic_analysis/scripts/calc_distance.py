import argparse
import datetime
import sys 
from collections import Counter
import itertools
import numpy as np
import pandas as pd
from collections import defaultdict


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def format_seqs_arr(s, n_seqs, nucs=np.array([97, 99, 103, 117])):
    n_seqs = int(n_seqs)
    size = int(len(s)/n_seqs)
    seqs_arr = \
        np.frombuffer(s.lower().encode(), dtype=np.int8)
    seqs_arr = seqs_arr.copy()
    if type(nucs) != type(None):
        seqs_arr[~np.in1d(seqs_arr, nucs)] = 110
    seqs_arr = \
        seqs_arr.reshape((n_seqs, int(seqs_arr.shape[0]/n_seqs)))
    return(seqs_arr)


def import_fasta(fasta_path, 
  nucs=np.array([97, 99, 103, 116])):
    s_names = []
    all_s = ''
    fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            all_s += s
    fh.close()
    s_arr = format_seqs_arr(all_s, len(s_names), nucs=nucs)
    s_names = np.array(s_names)
    return(s_arr, s_names)


def get_lineage_seqs(target_seq, ref_seq, exog_seqs):
    # any non ACTG- characters are converted to N
    # and Ns are not considered meaningful differences
    #exog_seqs_path = 'data/gisaid_hcov-19_2021_01_26_aligned_ref_filtered_masked.fasta'
    #target_seqs_path = 'data/Early_GA_43seqs_021121_wuhan_aligned_ref_filtered_masked_noref.fasta'
    # allows for ambiguos nucleotides
    target_poly_sites = \
        np.where(((target_seq != ref_seq) & 
            (target_seq != 110) & 
            (ref_seq != 110)))[0]
    target_poly_nucs = target_seq[target_poly_sites]
    # seqs in the same lineage but descendant (or match) of target seq
    # will match the target seq at all polymorphic positions
    # EXLCUDES sequences with Ns at these sites
    exog_descendant_arr_idxs = \
        np.where((exog_seqs[:,target_poly_sites] == target_poly_nucs).all(axis=1)==True)[0]
    # seqs in the same lineage but ancestral to target seq
    # will either match the reference or target 
    # at target poly sites and must match reference
    # or be N at all other sites
    # excludes sequences which are N at target_poly_sites
    mask = np.ones(ref_seq.shape[0], bool)
    mask[target_poly_sites] = False
    exog_ancestral_arr_idxs = \
        np.where(((((exog_seqs[:,mask]==ref_seq[mask]) | (exog_seqs[:,mask] == 110)).all(axis=1))
            & (((exog_seqs[:,target_poly_sites]==ref_seq[target_poly_sites]) | 
                (exog_seqs[:,target_poly_sites]==target_poly_nucs)).all(axis=1))) == True)[0]
    target_poly_nucs_str = list(target_poly_nucs.tobytes().decode("utf-8").upper())
    target_poly_sites_str = [str(i) for i in target_poly_sites]
    print(f'the query sequence has the following substitutions relative to the reference')
    print(list(zip(target_poly_sites_str, target_poly_nucs_str)))
    print(f'{(exog_seqs[:,target_poly_sites] != 110).all(axis=1).sum()} of the target region sequences were successfully genotyped at these positions')
    # need to identify exact matches and remove those 
    exog_exact_arr_idxs = \
        exog_exact_arr_idxs = \
        np.where((exog_seqs[:,target_poly_sites] == target_poly_nucs).all(axis=1) & 
            ((exog_seqs[:,mask] == target_seq[mask]) | 
                (exog_seqs[:,mask] == ref_seq[mask]) | 
                (exog_seqs[:,mask] == 110)).all(axis=1))[0]
    print(f'{exog_ancestral_arr_idxs.shape[0] - exog_exact_arr_idxs.shape[0]} of the target region sequences have a subset of the query sequence substitutions')
    print(f'they match the reference or are N at all other sites')
    print(f'{exog_exact_arr_idxs.shape[0]} of the target region sequences have all the query sequence substitutions')
    print(f'they match the reference or are N at all others')
    print(f'{exog_descendant_arr_idxs.shape[0] - exog_exact_arr_idxs.shape[0]} of of the target region sequences share the substitutions in the query sequence')
    print('they may have additional substitutions')
    lineage_seqs = np.unique(np.hstack([exog_ancestral_arr_idxs, exog_descendant_arr_idxs]))
    print(f'there are a total of {lineage_seqs.shape[0]} sequences in the lineage with the target sequence')
    return(lineage_seqs)


def run():
    parser = argparse.ArgumentParser()
    # input files
    # //// FASTA FILE PARSING ////
    parser.add_argument('--refSeq', default='EPI_ISL_402125.fasta')
    parser.add_argument('--seqs')
    parser.add_argument('--seqNameSep', default='|')
    parser.add_argument('--seqNameField', default=1)
    # //// METADATA FILE PARSING
    parser.add_argument('--metadata')
    parser.add_argument('--metadataDelim',
        default='\t')
    parser.add_argument('--metadataRegionCol',
        default=7, type=int)
    parser.add_argument('--metadataIDCol',
        default=1, type=int)
    parser.add_argument('--metadataDateCol',
        default=2, type=int)
    # //// SEQUENCE QUERY ////
    parser.add_argument('--querySeqName')
    parser.add_argument('--targetRegion')
    parser.add_argument('--minTargetSeqs',
        default=25, type=int)
    args = parser.parse_args()
    #args.seqs = 'for_steph/gisaid_hcov-19_2020_03_31_complete_hc_date_EHC_GA-EHC-069Q_metadata_aligned_ref_filtered_masked_imputed.fasta'
    #args.metadata = 'for_steph/metadata_aligned_GA-EHC-069Q.tsv'
    #args.querySeqName = '|GA-EHC-016P'
    #args.querySeqName = "|GA-EHC-031E"
    #args.targetRegion = 'MichiganUSA'
    #args.targetRegion = 'ColoradoUSA'
    #args.refSeq = 'for_steph/'+args.refSeq
    # import array not allowing for ambiguous nucleotides
    # ACTG only
    seqs_arr, seqs_names = import_fasta(args.seqs)
    ref_arr, ref_name = import_fasta(args.refSeq)
    ref_arr = ref_arr[0]
    seqs_names_idx_dict = {i:idx for idx, i in enumerate(seqs_names)}
    seq_name_dict = {i.split(args.seqNameSep)[args.seqNameField]: i for i in seqs_names}
    metadata = pd.read_csv(args.metadata, sep=args.metadataDelim, header=None)
    loc_dict = {seq_name_dict[i[args.metadataIDCol]]: i[args.metadataRegionCol] for i in metadata.values}
    date_dict = {seq_name_dict[i[args.metadataIDCol]]: i[args.metadataDateCol] for i in metadata.values}
    query_seq = seqs_arr[np.where(seqs_names==args.querySeqName)[0][0],:]
    if args.targetRegion:
        print(f'querying sequences from {args.targetRegion}')
        region_seqs_idxs = \
            metadata.loc[metadata[args.metadataRegionCol] == 
                args.targetRegion,args.metadataIDCol].\
                    map(seq_name_dict).\
                    map(seqs_names_idx_dict).values
    else:
        print(f'querying all input sequences')
        region_seqs_idxs = np.arange(seqs_arr.shape[0])
    region_seqs = \
        seqs_arr[region_seqs_idxs,:]
    region_seqs_names = \
        seqs_names[region_seqs_idxs]
    if len(region_seqs) < args.minTargetSeqs:
        print(f'less than {args.minTargetSeqs} sequences from {args.targetRegion}, results may be unreliable')
    print(f'query seq: {args.querySeqName}')
    print(f'target region: {args.targetRegion}')
    lineage_seqs_idxs = \
        get_lineage_seqs(query_seq, ref_arr, region_seqs)
    lineage_names = region_seqs_names[lineage_seqs_idxs]
    lineage_names = np.hstack([lineage_names, args.querySeqName])
    lineage_seqs = region_seqs[lineage_seqs_idxs]
    lineage_seqs = np.vstack([lineage_seqs, query_seq])
    lineage_subs = \
        np.where(((lineage_seqs != ref_arr) & 
            (ref_arr != 110) & 
            (lineage_seqs != 110)) == True)
    lineage_dists = np.bincount(lineage_subs[0])
    lineage_sub_sites = \
        np.split(lineage_subs[1], 
            np.cumsum(lineage_dists))[:-1]
    lineage_sub_ref_nucs = \
        [list(ref_arr[i].tobytes().decode('utf8').upper()) 
            for i in lineage_sub_sites]
    lineage_sub_nucs = \
        [list(lineage_seqs[idx, i].tobytes().decode('utf8').upper()) 
            for idx, i in enumerate(lineage_sub_sites)]
    lineage_subs_format = \
        [str(list(zip(lineage_sub_ref_nucs[idx], 
            lineage_sub_sites[idx]+1, 
            lineage_sub_nucs[idx]))) for idx, i in enumerate(lineage_sub_sites)]
    lineage_dates = [date_dict[i] for i in lineage_names]
    lineage_locs = [loc_dict[i] for i in lineage_names]
    out_df = pd.DataFrame([lineage_names, lineage_locs, 
        lineage_dates, lineage_dists, lineage_subs_format]).transpose().drop_duplicates()
    out_df = out_df.sort_values(3)
    out_df.to_csv(f'{args.querySeqName}_{args.targetRegion}_dists.tsv', sep='\t', header=None, index=None)




if __name__ == "__main__":
    run()



