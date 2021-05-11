import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import datetime
from utils import convert_to_datetime
from treetime.utils import numeric_date
import sys 
from collections import Counter
import itertools
from Bio.Seq import Seq


def format_seqs_arr(seqs_list):
    #seqs_arr = \
    #    np.fromstring(''.join([str(i.seq.lower()) for i in seqs_list]), dtype=np.int8)
    seqs_arr = \
        np.frombuffer(''.join([str(i.seq.lower()) for i in seqs_list]).encode(), dtype=np.int8)
    seqs_arr = np.copy(seqs_arr)
    # a = 97
    # c = 99
    # t = 116
    # g = 103
    # n = 110
    seqs_arr[(seqs_arr != 97) & 
        (seqs_arr != 99) & (seqs_arr != 116) & 
        (seqs_arr != 103)] = 110
    seqs_arr = \
        seqs_arr.reshape((len(seqs_list), len(seqs_list[0].seq)))
    return(seqs_arr)
    


def get_min_diffs(target_seqs_path, exog_seqs_path):
    #exog_seqs_path = 'non_ga_seqs.fasta'
    #target_seqs_path = 'ga_seqs.fasta'
    exog_seqs = list(SeqIO.parse(exog_seqs_path, 'fasta'))
    target_seqs = list(SeqIO.parse(target_seqs_path, 'fasta'))
    exog_seqs_arr = format_seqs_arr(exog_seqs)
    target_seqs_arr = format_seqs_arr(target_seqs)
    exog_target_dists = {}
    for exog_idx, exog_seq in enumerate(exog_seqs_arr):
        exog_target_dist = np.inf
        for target_seq_idx, target_seq in enumerate(target_seqs_arr):
            dist = ((exog_seq != target_seq) & 
                        (exog_seq != 110) & 
                        (target_seq != 110)).sum()
            exog_target_dist = \
                min(exog_target_dist, 
                    dist)
        exog_target_dists[exog_idx] = exog_target_dist
    # todo take name formatting as argument
    output = \
        pd.DataFrame([idx for idx,i in enumerate(exog_seqs)])
    output[1] = \
        output.loc[:,0].apply(lambda k: exog_seqs[k].description.split('|')[1])
    output[2] = output.loc[:,0].map(exog_target_dists)
    output = output[[1,2]]
    out_dist_path = \
        exog_seqs_path.replace('.fasta', '') + '_min_target_dist.tsv'
    output.to_csv(out_dist_path, header=None, index=None, sep='\t')
    weights = output[:]
    # using 1+ in case the distance is ever 0
    weight_vals = 1/(1+weights[2])
    weights[2] = 1/(1+weights[2])
    out_weight_path = \
        exog_seqs_path.replace('.fasta', '') + '_min_dist_weights.tsv'
    weights.to_csv(out_weight_path, header=None, index=None, sep='\t')


def run():
    parser = argparse.ArgumentParser()
    # input files
    # //// FASTA FILE PARSING ////
    parser.add_argument('--targetSeqs', 
        default=None, required=True)
    parser.add_argument('--exogSeqs', 
        default=None, required=True)
    args = parser.parse_args()
    get_min_diffs(args.targetSeqs, args.exogSeqs)



if __name__ == "__main__":
    run()



