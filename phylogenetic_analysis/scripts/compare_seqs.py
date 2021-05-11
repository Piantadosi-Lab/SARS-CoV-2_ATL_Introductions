import argparse
import pandas as pd
from Bio import SeqIO



def get_snp_aln(seqs=None, ref=None, 
  nuc_dict={'A':'A', 'C':'C', 'T':'T', 'G':'G', 'N': ['A', 'C','T','G']}, 
  out_name=None):
    # we want all the sites where at least one seq differs from wuhan
    snp_sites = set([])
    for seq in seqs:
        sequence = seq.seq.upper()
        remove_chars = '|'.join(set(sequence) - set(nuc_dict['N']) - set('N'))
        if len(remove_chars) > 0:
            sequence = re.sub(remove_chars, 'N', str(sequence))
        snp_sites.update([idx for idx in range(len(sequence)) if 
            (sequence[idx] != ref.seq[idx]) and 
            (sequence[idx] not in nuc_dict[ref.seq[idx]]) and 
            (ref.seq[idx] not in nuc_dict[sequence[idx]])])
    snp_sites = sorted(list(snp_sites))
    # gets nucleotide at each segregating site for each focal seq
    snp_alignment = \
        [[seq[i] for i in snp_sites] for seq in seqs]
    snp_alignment.insert(0,[ref.seq[i] for i in snp_sites])
    # convert to dataframe
    snp_alignment = \
        pd.DataFrame(snp_alignment, columns=snp_sites, index=[ref.description, *[i.description for i in seqs]])
    # make snp positions 1-indexed
    snp_alignment.columns = snp_alignment.columns + 1
    # want to get the distance matrix 
    print(out_name)
    snp_alignment.to_csv(out_name, sep='\t')
    return(snp_alignment)


def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--seqs', 
                        help='sequence alignment',
                        default=None)
    parser.add_argument('--seqNameSep', 
                        help='string to seperate sequence names',
                        default='|')
    parser.add_argument('--seqNameField', 
                        help='field in seperated sequence name',
                        default=1, type=int)
    parser.add_argument('--focalSeqs', 
                        help='list of focal sequences to focus on',
                        default=None)
    parser.add_argument('--focalSeqsDelim', 
                        default='\t')
    parser.add_argument('--focalSeqsCol', 
                        default=0, type=int)
    parser.add_argument('--focalSeqsSep', 
                        default='|')
    parser.add_argument('--focalSeqsField', 
                        default=1, type=int)
    parser.add_argument('--refSeq', 
                        help='reference fasta file',
                        default=None)
    args = parser.parse_args()


    nuc_dict = {"A": "A", "C": "C", "T": ["T", "U"], 
            "G": "G", "U": ["T", "U"], "R": ["A", "G"], 
            "Y": ["C", "T", "U"], "S": ["G", "C"], 
            "W": ["A", "T", "U"], "K": ["G", "T", "U"], 
            "M": ["A", "C"], "B": ["C", "G", "T", "U"], 
            "D": ["A", "G", "T", "U"], "H": ["A", "C", "T", "U"], 
            "V": ["A", "C", "G"], 
            "N": ["N", "A", "C", "T", "G", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]}
    seqs = list(SeqIO.parse(args.seqs, 'fasta'))
    ref = SeqIO.read(args.refSeq, 'fasta')
    focal_seq_names = \
        set(pd.read_csv(args.focalSeqs, header=None, 
            sep=args.focalSeqsDelim)[args.focalSeqsCol].str.split(args.focalSeqsSep, 
                expand=True)[args.focalSeqsField])
    focal_seqs = [i for i in seqs if 
        i.description.split(args.seqNameSep)[args.seqNameField] in focal_seq_names]
    print(f'{len(focal_seqs)} focal seqs found in alignment')
    # RETURNS 1 INDEXED COORIDNATES
    focal_snp_aln = \
        get_snp_aln(seqs=focal_seqs, ref=ref, nuc_dict=nuc_dict, 
            out_name='.'.join(args.focalSeqs.split('.')[:-1])+'_snp_aln.tsv')




if __name__ == "__main__":
    run()


