import numpy as np
import argparse


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


def get_maf(seqs, counted_nucs):
    from collections import Counter
    mask = np.zeros(seqs.shape, dtype=bool)
    bad_nucs = ~np.in1d(seqs, counted_nucs).reshape(seqs.shape)
    mask[bad_nucs] = True
    seqs[mask] = 110
    seqs_counted = [Counter(i) for i in seqs.T]
    # remove monomorphic sites, remove N
    seq_counted_filtered = \
        [i if len(i.keys()) > 1 else {} for 
            i in seqs_counted]
    _ = [i.pop(110, None) for 
        i in seq_counted_filtered]
    # sum all but largest value for each
    mac = \
        np.array([sum(sorted(list(i.values()))[:-1]) for 
            i in seq_counted_filtered])
    maf = mac/seqs.shape[0]
    return(np.array(maf))




def impute_nucs(seqs, target_seq_index, impute_sites, 
  impute_allow_nucs=[97, 99, 103, 116]):
    # get distance between target seq and all other seq
    target_seq = seqs[target_seq_index,:]
    dists = ((seqs != target_seq) & 
    	(target_seq != 110) & 
    	(seqs != 110)).sum(axis=1)
    # converg self-self comparison to inf so we can take min
    #dists[target_seq] = 1000000000
    # sort s_arr by distance
    sort_indices = np.argsort(dists)
    dists = dists[sort_indices]
    seqs = seqs[sort_indices, :]
    # get closest match at each nuceltoide
    impute_nucs = \
    	[[i for i in seqs[:,item] if 
    		i in impute_allow_nucs][0] for item in 
    		impute_sites]
    return(impute_nucs)



def write_seqs(chunked_s, names, path):
	with open(path, 'w') as fp:
		for idx,item in enumerate(chunked_s):
			split_seq = [item[i:i+80] for i in range(0, len(item), 80)]
			fp.write('>'+names[idx]+'\n')
			for chunk in split_seq:
				fp.write(chunk+'\n')


def export_seq_arr(seqs, seqs_names, out_path):
    raw_s = seqs.tobytes().decode("utf-8").upper()
    chunked_s = [raw_s[i:i+seqs.shape[1]] for i in range(0, len(raw_s), seqs.shape[1])]
    write_seqs(chunked_s, seqs_names, out_path)



def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seqs')
    parser.add_argument('--seqNameSep', default='|')
    parser.add_argument('--seqNameField', default=1)
    parser.add_argument('--ref', default='data/EPI_ISL_402125.fasta')
    parser.add_argument('--targetSeq')
    parser.add_argument('--mafThreshold', type=float, default=0.05)
    parser.add_argument('--outName', default='imputed')
    args = parser.parse_args()
    #args.seqs = 'for_steph/gisaid_hcov-19_2020_03_31_complete_hc_date_EHC_GA-EHC-069Q_metadata_aligned_ref_filtered_masked.fasta'
    #args.targetSeq = 'GA-EHC-069Q'
    num_nuc_dict = \
        {97: 'A', 110: 'N', 
            99: 'C', 116: 'T', 103: 'G'}
    s_arr, s_names = import_fasta(args.seqs, nucs=None)
    ref_arr, ref_name = import_fasta(args.ref, nucs=None)
    # get target seq
    # get distance between target seq and all other seqs
    s_names_format = \
        [i.split(args.seqNameSep)[args.seqNameField] for 
            i in s_names]
    target_seq_index = s_names_format.index(args.targetSeq)
    target_seq = s_arr[target_seq_index,:]
    # get MAF
    target_n_sites = np.where(target_seq == 110)[0]
    mafs = get_maf(s_arr[:,target_n_sites], 
        [np.frombuffer(i.encode(), dtype=np.int8) for i in ['actg']])
    # get sites with MAFs > threshold
    impute_sites = \
        target_n_sites[np.where(mafs > args.mafThreshold)]
    imputed_nucs = \
        impute_nucs(s_arr, target_seq_index, impute_sites, 
            impute_allow_nucs=np.frombuffer('actg'.encode(), dtype=np.int8))
    alts = np.array(imputed_nucs) == ref_arr[0,impute_sites]
    alts = ['Ref' if i==True else 'Alt' for i in alts]
    print('imputing the following nucleotides based on the closest related sequence')
    for idx,item in enumerate(imputed_nucs):
        print(f"at position {impute_sites[idx]+1} imputing nucleotide {num_nuc_dict[item]} ({alts[idx]})")
    s_arr[target_seq_index,impute_sites] = imputed_nucs
    export_seq_arr(s_arr, s_names, args.outName + '.fasta') 



if __name__ == '__main__':
    run()

