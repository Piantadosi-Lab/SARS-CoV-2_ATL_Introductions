import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re 
from collections import Counter

def get_shared_mutations(focal_seqs=None, ref=None, allow_chars=set(['A', 'C', 'T', 'G', 'N'])):
    n_focal_seqs = len(focal_seqs)
    snps = []
    ref_seq = list(ref.seq)
    for seq in focal_seqs:
        seq_seq = list(seq.seq)
        snps.extend([(seq.description, idx, ref[idx], i) for idx,i in \
            enumerate(seq_seq) if i != ref[idx] and i in allow_chars])
    snps = pd.DataFrame(snps, columns=['name', 'pos', 'ref', 'alt'])
    return(snps)


def calc_seg_sites(seqs, nuc_dict):
    # creats list of list and removes non-ACTG items
    seqs_lol = []
    for seq in seqs:
        sequence = seq.seq.upper()
        remove_chars = '|'.join(set(sequence) - set(nuc_dict['N']) - set('N'))
        if len(remove_chars) > 0:
            sequence = re.sub(remove_chars, 'N', str(sequence))
        seqs_lol.append(list(sequence))
    # transposes lol and gets set of nucleotides at each site
    seqs_lol_t = [set(i) for i in zip(*seqs_lol)]
    seg_sites = [[set(nuc_dict[k]) for k in i] for i in seqs_lol_t]
    seg_sites = [set.intersection(*i) for i in seg_sites]
    # segregating sites are those wiht no overlap in the sets
    seg_sites = [idx for idx,i in enumerate(seg_sites) if len(i) == 0]
    return(seg_sites)


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
    # convert to dataframe
    snp_alignment = \
        pd.DataFrame(snp_alignment, columns=snp_sites, index=[i.description for i in seqs])
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
    parser.add_argument('--metadata', 
                        help='metadata file',
                        default=None)
    parser.add_argument('--metadataDelim', 
                        help='metadata file delimiter',
                        default='\t')
    parser.add_argument('--metadataIDCol', 
                        help='metadata ID column',
                        default=1, type=int)
    parser.add_argument('--metadataDateCol', 
                        help='metadata ID column',
                        default=2, type=int)
    parser.add_argument('--outName', 
                        help='outname for file',
                        default=None)
    args = parser.parse_args()

    #args.seqs = 'data/gisaid_hcov-19_2020_03_31_complete_hc_date_EHC_metadata_aligned_ref_filtered_masked.fasta'
    #args.focalSeqs = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0_family.tsv'
    #args.refSeq = 'data/EPI_ISL_402125.fasta'
    #args.outName = 'test'
    #args.metadata = 'data/metadata_aligned.tsv'

    nuc_dict = {"A": "A", "C": "C", "T": ["T", "U"], 
            "G": "G", "U": ["T", "U"], "R": ["A", "G"], 
            "Y": ["C", "T", "U"], "S": ["G", "C"], 
            "W": ["A", "T", "U"], "K": ["G", "T", "U"], 
            "M": ["A", "C"], "B": ["C", "G", "T", "U"], 
            "D": ["A", "G", "T", "U"], "H": ["A", "C", "T", "U"], 
            "V": ["A", "C", "G"], 
            "N": ["N", "A", "C", "T", "G", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]}
        # contains all "compmlete" "high coverage" "sampling date complete" genomes
        # on GISAID sampled thorugh January 26 2021 with the L84S substitution
        # which is a clade defining SNP of 19B
        # plus Wuhan-Hu-1
    #args.seqs = 'data/weighted_downsampling_old/19B/l84s_EHC_masked_aligned.fasta'
    #args.focalSeqs = 'data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3_0.tsv'
    #args.ref = 'EPI_ISL_402125'
    # get focal seqs
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
            out_name=args.outName.replace('.tsv','')+'_snp_aln.tsv')
    # we want to get the shared SNPs, thus columns with only 1 value
    alleles = focal_snp_aln.apply(set, axis=0).tolist()
    # remove any character not in nuc_dict kids
    alleles = [list(i) for i in alleles]
    _ = [[i.remove(k) for k in i if k not in nuc_dict.keys()] 
            for i in alleles]
    alleles = [set(i) for i in alleles]
    # gets non-segregating sites
    shared_snps = \
        [focal_snp_aln.columns[idx] for idx, i in enumerate(alleles) if 
            len(set.intersection(*[set(nuc_dict[k]) for k in i])) > 0]
    # shared SNPs have at least 1 overlap
    print(f'there are {len(shared_snps)} shared SNPs between the {len(focal_seqs)} \
focal sequences: {" ".join([str(i) for i in shared_snps])}')
    mutation_profile = focal_snp_aln.loc[:,shared_snps].apply(
        lambda k: list(Counter(k).keys())[0]).to_frame()
    mutation_profile[1] = \
        pd.Series([ref.seq[i-1] for i in mutation_profile.index],
            index = mutation_profile.index)
    mutation_profile = mutation_profile[[1,0]]
    vcf_str = ['##fileformat=VCFv4.0\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n']
    vcf_str.extend([['NC_045512.2\t'+str(idx)+'\t.\t'+
            i[1]+'\t'+i[0]+'\t.\t.\t.\n'] for 
        idx, i in mutation_profile.iterrows()])
    if args.outName:
        mutation_profile.to_csv(args.outName+'_focal_snps_shared.tsv', 
            sep='\t', header=None)
        with open(args.outName+'_focal_snps_shared.vcf', 'w') as fp:
            for line in vcf_str:
                fp.writelines(line)
    else:
        mutation_profile.to_csv(
            args.focalSeqs.replace('.tsv', '')+'_focal_snps_shared.tsv', 
                sep='\t', header=None)
        with open(args.focalSeqs.replace('.tsv', '')+'_focal_snps_shared.vcf', 'w') as fp:
            for line in vcf_str:
                fp.writelines(line)
    # now that we have the mutational profile we can get the sequences in the alignment that match 
    # adjusting index by 1 so they are 0 indexed
    mutation_profile_dict = \
        {idx-1: i[0] for idx, i in mutation_profile.iterrows()}
    match_profile_seqs = \
        [i for i in seqs if 
            all([i[k] == mutation_profile_dict[k] for 
                k in mutation_profile_dict.keys()])]
    match_profile = \
        [i.description for i in match_profile_seqs]
    match_profile_split_dict = {i.split(args.seqNameSep)[args.seqNameField]: i for i in match_profile}
    match_profile_df = pd.DataFrame(match_profile)
    if args.metadata:
        metadata = pd.read_csv(args.metadata, header=None, sep=args.metadataDelim)
        match_profile_metadata = metadata[metadata[args.metadataIDCol].isin(match_profile_split_dict.keys())]
        id_col = match_profile_metadata.loc[:,args.metadataIDCol].map(match_profile_split_dict).copy(deep=True)
        match_profile_metadata.loc[:,match_profile_metadata.shape[1]] = id_col
        match_profile_df = match_profile_metadata
        match_profile_df = match_profile_df[[match_profile_df.shape[1]-1, *range(match_profile_df.shape[1]-1)]]
        print('dt start')
        dt = pd.to_datetime(match_profile_df.loc[:,args.metadataDateCol]).copy(deep=True)
        match_profile_df.loc[:,args.metadataDateCol] = dt
        match_profile_df = match_profile_df.sort_values(by=args.metadataDateCol)
        print(match_profile_df)
    if args.outName:
        with open(args.outName+'.fasta', 'w') as fp:
            SeqIO.write(match_profile_seqs, fp, 'fasta')
        match_profile_df.to_csv(
            args.outName+'_match_focal_snps.tsv', sep='\t', header=None, index=None)
    else:
        with open(args.focalSeqs.replace('.tsv','')+'_match_focal_snps.fasta', 'w') as fp:
            SeqIO.write(seqs, fp, 'fasta')
        match_profile_df.to_csv(
            args.focalSeqs.replace('.tsv','')+'_match_focal_snps.tsv', sep='\t', header=None, index=None)



if __name__ == "__main__":
    run()


