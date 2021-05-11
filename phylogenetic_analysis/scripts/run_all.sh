
# global sequence set includes everything marked
# "complete," "high coverage," and "collection date complete"
# from human host 
# on GISAID as of 2021-03-27
# sampled through 2020-03-31
globalSeqs=data/gisaid_hcov-19_2020_03_31_complete_hc_date.fasta
# GISAID global metadta
# downloaded on 2021-03-27
globalData=data/metadata.tsv
# EHC sequences
ehcSeqs=data/Early_GA_seqs.fasta
# EHC metadata
ehcData=data/Early_GA_seqs_metadata.tsv
# TODO update this
l84sSeqs=data/gisaid_hcov-19_2020_03_31_l84s_complete_hc_date.fasta
# reference information
refSeq='data/EPI_ISL_402125.fasta'
refName='EPI_ISL_402125'

# sets up folder structure
mkdir data
mkdir config
mkdir data/weighted_downsampling
mkdir data/19B_subclade
mkdir data/l84s
mkdir data/temporal_downsampling

# //// 1. BASIC FORMATTING OF EHC SEQUENCES AND METADATA ////
# formats the ehc_seqs header so they can be split on a "|"
# also removes "amp" and "merged" from file names
sed '/\>/ s//\>\|/g' $ehcSeqs | \
    sed 's/\_amp//g' | \
    sed 's/\_merged//g' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' > ${ehcSeqs%.fasta}_format.fasta 
# gets names which appear in both ehcSeqs
# and in globalSeqs (ehcSeqs already uploaded to gisaid)
# and adds them to exclude folder
# so we don't double count them when generating our downsampled alignment
rm -rf config/exclude.tsv
grep '>' ${ehcSeqs%.fasta}_format.fasta | \
    sed 's/\>|//g'  | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' | \
    grep -f /dev/stdin <(grep '>' $globalSeqs) | \
    awk -F'|' '{print "1\t"$2}' > config/exclude.tsv
# the following samples are replicates from the same patient, so we exclude one of them
echo '# same patient samples
# preferentially including NP samples
1\tGA-EHC-087I
1\tGA-EHC-081C' \
>> config/exclude.tsv
# formats the EHC metadata
tail -n +2 $ehcData | \
        sed 's/EHC_C19_/GA-EHC-/g' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' |\
        awk -F"\t" '{print $1"\t"$1"\t"$5"\tNorth America\tUSA\tGeorgia\tGA\tGeorgiaUSA"}' \
        > ${ehcData%.tsv}_format.tsv

# subsets metadata to just the sequences we have in fasta file
awk -F'\t' 'FNR==NR {hash[$1]; next} $1 in hash' \
    <(grep ">" ${ehcSeqs%.fasta}_format.fasta |\
    sed 's/\>\|//g') \
    ${ehcData%.tsv}_format.tsv \
    > ${ehcData%.tsv}_format_filter.tsv


# //// 2. GENERATE A SINGLE FORMATTED METADATA FILE FOR BOTH EHC AND GISAID DATA
# recomment this section!!!
# REMOVE AMBIGUOUS DATES 
# remove "exlcude" sequences
# human samples only
# no travel history
# no ambiguous dates
# no cruise ships
# complete and high quality (only these ones are downloaded)

awk -F'\t' 'FNR==NR {hash[$2]; next} !($2 in hash)' \
    <(cat config/exclude.tsv) \
    <(awk -F'\t' '{if ($15 == "Human" || $15 == "human") print $0}' $globalData | \
        awk -F'\t' '{if ($6$7$8 == $10$11$12) print $0}' | \
        awk -F'\t' '{ if ($5 !~ "X" && $5 !~ "x") print $0}' | \
        awk -F'\t' '{split($5,a,"-"); if (length(a)==3) print $0}' | \
        grep -i -v cruise | \
        awk -F'\t' '{print $1"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8}' | \
        sed "s/'/\_/g") \
    > ${globalData%.tsv}_filter.tsv


cat \
<(awk -F'\t' '{if ($5!="USA") print $0"\tNON_GA\t"$5}' ${globalData%.tsv}_filter.tsv) \
<(awk -F'\t' '{if ($5=="USA" && $6!="Georgia") print $0"\tNON_GA\t"$6"USA"}' ${globalData%.tsv}_filter.tsv) \
<(awk -F'\t' '{if ($5=="USA" && $6=="Georgia") print $0"\tGA\t"$6"USA"}' ${globalData%.tsv}_filter.tsv) \
> ${globalData%.tsv}_filter_format.tsv


# //// 2. CONCATENATES GLOBAL SEQUENCES AND EHC SEQUENCES THEN ALIGNS ////
cat $globalSeqs ${ehcSeqs%.fasta}_format.fasta \
    > ${globalSeqs%.fasta}_EHC.fasta
cat ${globalData%.tsv}_filter_format.tsv ${ehcData%.tsv}_format_filter.tsv > data/metadata_combined.tsv
# makes sure we only include sequences in the metadata
python3 scripts/get_seqs.py \
    --seqs ${globalSeqs%.fasta}_EHC.fasta \
    --getSeqs data/metadata_combined.tsv \
    --getSeqsCol 1 \
    --getSeqsField 0 \
    --outName ${globalSeqs%.fasta}_EHC_metadata

python3 scripts/align_seqs.py \
    --sequences ${globalSeqs%.fasta}_EHC_metadata.fasta \
    --referenceName $refName \
    --maskHead 100 \
    --maskTail 100 \
    --maskSites 11083 15324 21575  \
    --alignType reference \
    --minLength 29000

# subset the metadata to just the aligned sequences
# and combine with lineage assignments
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_combined.tsv \
    > ${globalData%.tsv}_aligned.tsv


# //// 3. CALCULATE DOWNSAMPLING WEIGHTS ////
# downsampling scheme is two tiered: 
# 1) weightedd on a per country basis related to the number of 
# cases in each country by March 31st
# 2) 1/(D+1) where D is the maximum distance 
# between that sequence and all of the GA sequences
# // 3A. COUNTRY CASE WEIGHTS //
# TODO convert data to tsv for consistencies (doesn't affect results at all)
python3 scripts/calc_case_weights.py \
    --globalDat data/time_series_covid19_confirmed_global.csv \
    --usDat data/time_series_covid19_confirmed_US.csv \
    --dateMax 03/31/20 \
    --focalStates Georgia


# // 3B. DISTANCE WEIGHTS //
# first get all GA sequences from global alignment
# then concat that with EHC alignment
# NOTE: THIS FILE WILL DOUBLE COUNT THOSE EHC SEQUENCES
# WHICH ARE ALREADY ON GISAID. THIS IS FINE 
# FOR THE PURPOSES OF CALCULATING WEIGHTS BUT THE 
# GA SEQS FILE SHOULD NOT BE USED FOR ANY OTHER PURPOSES. 
# while metadata includes information for *everything* on GISAID
# the alignment only contains sequences through march 31
# and which pass filtering, so "gisaid_ga_seq_ids.txt"
# will include sequences which will not end up in teh 
# ga_seqs.fasta file
grep 'USA\tGeorgia' ${globalData%.tsv}_aligned.tsv | \
    awk '{print $2 }' > data/ga_seq_ids.txt
# get GA seqs
python3 scripts/get_seqs.py \
    --outName data/ga_seqs \
    --seqs ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --getSeqs data/ga_seq_ids.txt \
    --getSeqsField 0
# calculate minimum distance
python3 scripts/calc_min_distance.py \
    --targetSeqs data/ga_seqs.fasta \
    --exogSeqs ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta

# //// 4. DOWNSAMPLE THE GLOBAL SEQUENCES BASED ON MINIMUM DISTANCE ////
python3 scripts/downsample_seqs.py \
    --sequences ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --include ./config/include.tsv \
    --exclude ./config/exclude.tsv \
    --maxDate 2020-03-31 \
    --minDate 2019-01-01 \
    --regionWeights data/country_case_weights.tsv \
    --interRegionWeights ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked_min_dist_weights.tsv \
    --targetN 6000 \
    --outName data/weighted_downsampling/ga_focused_aligned_masked_weighted

# get clade assignments for downsampled sequences
nextclade \
    --input-fasta data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --output-tsv data/weighted_downsampling/ga_focused_aligned_masked_weighted_clades.tsv

# get metadata for all included sequences and combine with clades
awk -F'\t' 'FNR==NR {split($1, b, "|"); a[b[2]]=$2; next}{print $0"\t"a[$2]}' \
    data/weighted_downsampling/ga_focused_aligned_masked_weighted_clades.tsv \
    <(awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
        <(grep ">" data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta |\
        sed 's/\>//g' | awk -F'|' '{print $2}') \
        ${globalData%.tsv}_aligned.tsv) \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv
# get metadata for just GA included sequences
grep "\tGA\t" data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv


# //// 5. BUILD ML TREE OF WEIGHTED DOWNSAMPLED SEQS ////
iqtree2 -redo --polytomy -m TEST -T AUTO -bb 1000 --wbtl -T 4 --prefix data/weighted_downsampling/ga_focused_aligned_masked_weighted -s data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta
iqtree -redo --polytomy -m TEST -T AUTO -bb 1000 --wbtl -T 4 --prefix ga_focused_aligned_masked_weighted -s ga_focused_aligned_masked_weighted.fasta


# //// 6. ESTIMATE NUMBER OF IMPORTATIONS BASED ON WEIGHTED ML TREE ////
python3 scripts/estimate_importations.py \
    --trees data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile \
    --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    --sequences data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --outgroup $refName \
    --biasCorrection 2.5 \
    --regions GA \
    --metadataLocCol 6 \
    --saveAln 

# plot RTT
python3 scripts/plot_rtt.py \
    --data data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_rtt.csv \
    --outName figures/weighted_downsampling_rtt

# -------------------------------------------------------------------#
# //// 7. REPEAT ON 100 BOOTSTRAP REPLICATES OF WEIGHTED ML TREE ////
refName='EPI_ISL_402125'
python3 scripts/estimate_importations.py \
            --trees data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot \
            --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
            --sequences data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
            --outgroup $refName \
            --biasCorrection 2.5 \
            --regions GA \
            --metadataLocCol 6 \
            --bootstrapReps 0 100

# -------------------------------------------------------------------#
# //// 8. RESULTS FROM STEPS 7 AND 8 ////
python3 scripts/plot_annotated_tree.py \
    --tree data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick \
    --treeStates data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_node_states.csv \
    --treeBootstrapDir "data/weighted_downsampling_old/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv" \
    --config config/annotated_tree.json \
    --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    --travelData data/ga_travel_seqs.tsv
    # get # of introductions
python3 scripts/estimate_introductions_ml.py \
    --treeBootstrapDir "data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv"

# plot ML divergence tree
python3 scripts/plot_divergence_tree.py \
    --tree data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined.newick\
    --config config/annotated_global_divergence_tree.json \
    --metadata data/metadata_combined.tsv \
    --metadataRegionCol 6


# --------------------------------------------------------------------------------------#
# //// 8. GET MUTATIONAL PROFILE SHARED BETWEEN GRANDPARENT OF MRCA OF 19B SUBCLADE ////
# cluster GA seqs in tree
python3 scripts/cluster_seqs.py \
    --tree  data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick  \
    --seqNames data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv 
# subset to just cluster 0
# always the biggest cluster, they're sorted
grep '\t0' data/weighted_downsampling/ga_focused_aligned_masked_weighted_clusters_0.3.tsv | \
    awk '{print $1}' \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv
# get the mutational profile of these sequences
python3 scripts/match_mutational_profile.py \
    --seqs  data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --focalSeqs data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv \
    --refSeq $refSeq \
    --outName data/weighted_downsampling/19B_subclade
# annotate the vcf
java -jar scripts/snpEff/snpEff.jar \
    -noInteraction -no-upstream -no-downstream -no-intergenic \
    NC_045512.2 data/weighted_downsampling/19B_subclade_focal_snps_shared.vcf \
    > data/weighted_downsampling/19B_subclade_focal_snps_shared_ann.vcf
# get all sequence names in the tree which descend from the grandparent of the MRCA of these sequences
python3 scripts/get_seq_names_from_tree.py \
    --tipNames data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv \
    --tree data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick \
    --nGenerations 3


# --------------------------------------------------------------------------------------#
# //// 9. GET ALL AVAILABLE SEQUENCES WHICH MATCH THIS MUTATIONAL PROFILE ////
# --------------------------------------------------------------------------------------#
# get the mutational profile of all sequences in this cluster family
# and pull all sequences out of the global data matching this mutational profile
python3 scripts/match_mutational_profile.py \
    --seqs  ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --focalSeqs data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0_family.tsv \
    --refSeq data/$refName.fasta \
    --outName data/19B_subclade/19B_subclade_family
# get metadata
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/19B_subclade/19B_subclade_family.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_aligned.tsv \
    > data/19B_subclade/19B_subclade_family.tsv
# build ML tree
iqtree2 -redo --polytomy -m TEST -T AUTO --wbtl -T 4 \
    --prefix data/19B_subclade/19B_subclade_family \
    -s data/19B_subclade/19B_subclade_family.fasta
# run this through treetimes clock filter
# also reroots the tree using least squares 
# and outputs data for RTT plot
python3 scripts/clock_filter.py \
    --tree data/19B_subclade/19B_subclade_family.treefile \
    --seqs data/19B_subclade/19B_subclade_family.fasta \
    --metadata data/19B_subclade/19B_subclade_family.tsv \
    --metadataDateCol 2 \
    --metadataIDCol 1 \
    --nIQD 4
# get metadata
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/19B_subclade/19B_subclade_family_clockfilter.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_aligned.tsv \
    > data/19B_subclade/19B_subclade_family_clockfilter.tsv
# rtt plot
python3 scripts/plot_rtt.py \
    --data data/19B_subclade/19B_subclade_family_rtt.csv\
    --outName figures/19B_subclade_rtt
# ml plot
python3 scripts/plot_divergence_tree.py \
    --tree data/19B_subclade/19B_subclade_family_clockfilter.newick \
    --config config/annotated_divergence_tree.json \
    --metadata data/19B_subclade/19B_subclade_family.tsv  \
    --metadataRegionCol 7

# --------------------------------------------------------------------------------------#
# //// 9. GENERATE BEAST XML FILE AND RUN////
# subset the metadata to just the samples we want
# 1. U.S. states with more than 3 samples
# Exclude sequences with no USA state data
# 2. international sequences sampled before March
# there are no 2019 samples in this clade, so just filter by month
awk -F '\t' '{if ($5=="USA") print $6}' data/19B_subclade/19B_subclade_family_clockfilter.tsv | sort | \
    uniq -c | awk '{if ($1 > 3) print $2 " " $3}' |  grep -v USA | sed -e 's/[[:space:]]*$//' | \
    sed '/^$/d' > data/19B_subclade/include_states.tsv

cat <(awk -F'\t' '{if ($5=="USA") print $0}' data/19B_subclade/19B_subclade_family_clockfilter.tsv | \
        grep -f data/19B_subclade/include_states.tsv /dev/stdin) \
    <(awk -F'\t' '{if ($5 != "USA") print $0}' data/19B_subclade/19B_subclade_family_clockfilter.tsv | \
        awk -F'\t' '{split($3,a,"-"); if (a[2] < 3) print $0}') \
    > data/19B_subclade/19B_subclade_beast_include.tsv


python3 scripts/generate_xml.py \
    --xmlTemplate config/beast2_trait_template.xml \
    --outName data/19B_subclade/19B_subclade \
    --seqs data/19B_subclade/19B_subclade_family_clockfilter.fasta \
    --metadata data/19B_subclade/19B_subclade_beast_include.tsv \
    --metadataTraitCol 7 \
    --alnName 19B_sublcade \
    --traitName location 



# run beast
# generate MCC tree
/Applications/BEAST\ 2.6.3/bin/treeannotator -heights median -b 10  -lowMem \
    data/19B_subclade/19B_sublcade_location_tree_with_trait.trees \
    data/19B_subclade/19B_location_mcc.tre


# parse set of beast trees
python3 scripts/parse_beast_trees.py \
    --trees data/19B_subclade/19B_sublcade_location_tree_with_trait.trees 

# gets statistics
# todo make this more efficient!
python3 scripts/estimate_importations_beast.py \
    --trees data/19B_subclade/19B_sublcade_location_tree_with_trait.pkl \
    --focalRegion GeorgiaUSA


# --------------------------------------------------------------------------------------#
# //// 11. GET FUTURE SEQUENCES WHICH MATCH MUTATIONAL PROFILE of GA 19B subclade ////
# need to align and mask L84S sequences
# contains all "compmlete" "high coverage" "sampling date complete" genomes
# on GISAID sampled thorugh January 26 2021 with the L84S substitution
# which is a clade defining SNP of 19B
# add ref for aligning
cat $refSeq $l84sSeqs > ${l84sSeqs%.fasta}_wuhan.fasta
# get all 19B EHC seqs 
grep 'GA-EHC' data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv \
    > data/l84s/19B_subclade_ehc.tsv

python3 scripts/get_seqs.py \
    --seqs data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --getSeqs data/l84s/19B_subclade_ehc.tsv\
    --getSeqsCol 0 \
    --outName data/l84s/19B_sublcade_EHC
# add in 19B subclade EHC seqs
cat ${l84sSeqs%.fasta}_wuhan.fasta \
    data/l84s/19B_sublcade_EHC.fasta \
    > ${l84sSeqs%.fasta}_wuhan_EHC.fasta

# makes sure we only include sequences in the metadata
# also removes extra Emory sequences
python3 scripts/get_seqs.py \
    --seqs ${l84sSeqs%.fasta}_wuhan_EHC.fasta \
    --getSeqs data/metadata_combined.tsv \
    --getSeqsCol 1 \
    --getSeqsField 0 \
    --outName data/l84s/l84s_filtered

# align them
python3 scripts/align_seqs.py \
    --sequences data/l84s/l84s_filtered.fasta \
    --referenceName $refName \
    --maskHead 100 \
    --maskTail 100 \
    --maskSites 11083 15324 21575  \
    --alignType reference \
    --minLength 29000

# get the metadata
# subset the metadata to just the aligned sequences
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/l84s/l84s_filtered_aligned_ref_filtered_masked_noref.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_combined.tsv \
    > data/l84s/l84s_aligned.tsv

# FINALLY get sequences which match mutational profile
python3 scripts/match_mutational_profile.py \
    --seqs data/l84s/l84s_filtered_aligned_ref_filtered_masked_noref.fasta \
    --focalSeqs data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv \
    --refSeq $refSeq \
    --metadata data/l84s/l84s_aligned.tsv \
    --outName data/l84s/l84s_19B_subclade

# get two earliest sequences
head -n 2 data/l84s/l84s_19B_subclade_match_focal_snps.tsv | \
    awk -F'\t' '{print $3}'\
    > data/l84s/l84s_19B_subclade_earliest.tsv
python3 scripts/compare_seqs.py \
    --seqs data/l84s/l84s_filtered_aligned_ref_filtered_masked_noref.fasta \
    --focalSeqs data/l84s/l84s_19B_subclade_earliest.tsv \
    --refSeq data/EPI_ISL_402125.fasta \
    --focalSeqsField 0



# //// 12. GENERATE BEAST/SUBCLADE PLOTS ////

# beast plot
python3 scripts/plot_annotated_beast.py \
    --tree data/19B_subclade/19B_location_mcc.tre \
    --mutationalProfile data/19B_subclade/19B_subclade_EHC_focal_snps_shared.tsv \
    --config config/annotated_beast_tree.json \
    --matchMutationalProfile data/l84s/l84s_19B_subclade_match_focal_snps.tsv \
    --treeFocusTips data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv


# //// 13. DOWNSAMPLE THE GLOBAL SEQUENCES BASED ON SAMPLING DATE ////
python3 scripts/downsample_seqs.py \
    --sequences ${globalSeqs%.fasta}_EHC_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --include ./config/include.tsv \
    --exclude ./config/exclude.tsv \
    --maxDate 2020-03-31 \
    --samplesPerWeek 20 \
    --outName data/temporal_downsampling/ga_focused_aligned_masked_temporal

# subset metadata
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_aligned.tsv \
    > data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv
# get metadata for just GA included sequences
grep "\tGA\t" data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
    > data/temporal_downsampling/ga_focused_aligned_masked_temporal_ga_included_seqs.tsv

# //// 14. BUILD ML TREE OF TEMPORAL DOWNSAMPLING SEQUENCES ////
iqtree2 -redo --polytomy -m TEST -T AUTO -bb 1000 --wbtl -T 4 --prefix data/temporal_downsampling/ga_focused_aligned_masked_temporal -s data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta

# //// 15. ESTIMATE NUMBER OF IMPORTATIONS BASED ON TEMPORAL ML TREE ////
python3 scripts/estimate_importations.py \
            --trees data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile \
            --metadata data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
            --sequences data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta \
            --outgroup $refName \
            --biasCorrection 2.5 \
            --regions GA \
            --metadataLocCol 6 \
            --saveAln 

# //// 16. REPEAT ON 100 BOOTSTRAP REPLICATES OF TEMPORAL ML TREE ////
refName='EPI_ISL_402125'
python3 scripts/estimate_importations.py \
            --trees data/temporal_downsampling/ga_focused_aligned_masked_temporal.ufboot \
            --metadata data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
            --sequences data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta \
            --outgroup $refName \
            --biasCorrection 2.5 \
            --regions GA \
            --metadataLocCol 6 \
            --bootstrapReps 0 100

# //// 17. RESULTS FROM STEPS 15 AND 16 ////
python3 scripts/plot_annotated_tree.py \
    --tree data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_refined_time.newick \
    --treeStates data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_refined_node_states.csv \
    --treeBootstrapDir "data/temporal_downsampling/ga_focused_aligned_masked_temporal.ufboot_tres/*/*_importations.csv" \
    --config config/annotated_temporal_tree.json \
    --metadata data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
    --travelData data/ga_travel_seqs.tsv

# plot RTT
python3 scripts/plot_rtt.py \
    --data data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_rtt.csv \
    --outName figures/temporal_downsampling_rtt
    
# get # of introductions
python3 scripts/estimate_introductions_ml.py \
    --treeBootstrapDir "data/temporal_downsampling/ga_focused_aligned_masked_temporal.ufboot_tres/*/*_importations.csv" 


# //// GENERATE GISAID ACKNOWLEDGEMENT TABLE ////
# get all epi ISLs used
cat <(awk -F'\t' '{print $2}' data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv | grep -v "GA-EHC") \
    <(awk -F'\t' '{print $2}' data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv | grep -v "GA-EHC") \
    <(awk -F'\t' '{print $2}' data/l84s/l84s_aligned.tsv | grep -v "GA-EHC") \
    <(awk -F'\t' '{print $2}' config/exclude.tsv | grep -v "#"  | grep -v "GA-EHC" | grep -v "^$") | \
    uniq > \
    data/all_included_seqs.txt


# //// NEW STUFF ENDS HERE ////



# TEST RUN RECONSTRUCTING GA/EHC
python3 scripts/estimate_importations.py \
    --trees data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_test \
    --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs_test.tsv \
    --sequences data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --outgroup $refName \
    --biasCorrection 2.5 \
    --regions GeorgiaUSA GeorgiaEHCUSA \
    --metadataLocCol 7 \
    --saveAln 

python3 scripts/plot_annotated_tree.py \
    --tree data/weighted_downsampling/old/ga_focused_aligned_masked_weighted.treefile_test_tres/0/0_refined_time.newick \
    --treeStates data/weighted_downsampling/old/ga_focused_aligned_masked_weighted.treefile_test_tres/0/0_refined_node_states.csv \
    --treeBootstrapDir "data/weighted_downsampling/old/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv" \
    --config config/annotated_tree_test.json \
    --metadata data/weighted_downsampling/old/ga_focused_aligned_masked_weighted_all_included_seqs.tsv


# //// 8. GET 19B SEQUENCES WHICH PASSED CLOCK FILTER
# todo change from "old" to current
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_aligned.tsv \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_time.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv | \
    grep "19B" \
    > data/weighted_downsampling/19B/19B_seqs.tsv
python3 scripts/get_seqs.py \
    --outName data/weighted_downsampling/19B/19B \
    --seqs data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_time.fasta \
    --getSeqs data/weighted_downsampling/19B/19B_seqs.tsv \
    --getSeqsCol 1 \
    --getSeqsField 0

# //// 9. GENERATE BEAST XML FILE AND RUN////
# todo change from "old" to current
python3 scripts/generate_xml.py \
    --xmlTemplate config/beast2_trait_template.xml \
    --outName data/weighted_downsampling/19B/19B \
    --seqs data/weighted_downsampling/19B/19B.fasta \
    --metadata data/weighted_downsampling/19B/19B_seqs.tsv \
    --metadataTraitColPrimary 7 \
    --metadataTraitColSecondary 3 \
    --minTraitSize 1000000 \
    --alnName 19B \
    --traitName location
# run beast command
# treeannotator
/Applications/BEAST\ 2.6.3/bin/treeannotator -heights median -b 50  -lowMem \
    data/weighted_downsampling19B/location_tree_with_trait.trees \
    data/weighted_downsampling_old/19B/19B_location_mcc.tre

# //// 10. CLUSTER SEQUENCES IN BEAST MCC TREE
python3 scripts/cluster_seqs.py \
    --tree  data/weighted_downsampling_old/19B/19B_location_mcc.tre \
    --seqNames data/weighted_downsampling_old/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv \
# get sequence names in cluster 0 and format them
grep '\t0' data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3.tsv | \
    awk -F'\t' '{print $1}' | \
    sed 's/GA\_EHC\_/GA\-EHC\-/g' \
    > data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3_0.tsv

# //// 11. GET MUTATIONAL PROFILE OF GA 19B SEQS ////
# need to align and mask L84S sequences
# contains all "compmlete" "high coverage" "sampling date complete" genomes
# on GISAID sampled thorugh January 26 2021 with the L84S substitution
# which is a clade defining SNP of 19B
# add ref
cat data/$refName.fasta $l84sSeqs > ${l84sSeqs%.fasta}_wuhan.fasta
# subset to only things in metadata
# this removes non human samples and travel associated samples
python3 scripts/get_seqs.py \
    --seqs ${l84sSeqs%.fasta}_wuhan.fasta \
    --getSeqs ${globalData%.tsv}_combined.tsv \
    --getSeqsCol 1 \
    --outName ${l84sSeqs%.fasta}_wuhan_filtered
# explicitly remove any GA 19B sequences (we will readd later)
# this is to account fro the fact that not all EHC sequences are on GISAID
# also move over to 19B subfolder
python3 scripts/get_seqs.py \
    --seqs ${l84sSeqs%.fasta}_wuhan_filtered.fasta \
    --getSeqs data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3.tsv \
    --getSeqsCol 0 \
    --outName data/weighted_downsampling_old/19B/l84s_wuhan_filtered_noGA \
    --excludeSeqs
# align them
python3 scripts/align_seqs.py \
    --sequences data/weighted_downsampling_old/19B/l84s_wuhan_filtered_noGA.fasta \
    --referenceName $refName \
    --maskHead 100 \
    --maskTail 100 \
    --maskSites 11083 15324 21575  \
    --alignType reference \
    --minLength 29000
# readd already aligned GA 19B seqs
# need to format getSeqs file 
python3 scripts/get_seqs.py \
    --seqs data/weighted_downsampling_old/19B/19B.fasta \
    --getSeqs data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3.tsv \
    --getSeqsCol 0 \
    --outName data/weighted_downsampling_old/19B/19B_ga \
    --formatNames "GA_EHC_" "GA-EHC-"


cat data/weighted_downsampling_old/19B/19B_ga.fasta  \
    data/weighted_downsampling_old/19B/l84s_wuhan_filtered_noGA_aligned_ref_filtered_masked_noref.fasta\
    > data/weighted_downsampling_old/19B/l84s_EHC_masked_aligned.fasta

# FINALLY get mutational profile and sequences which match
python3 scripts/match_mutational_profile.py \
    --seqs data/weighted_downsampling_old/19B/l84s_EHC_masked_aligned.fasta \
    --focalSeqs data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3_0.tsv \
    --refSeq data/EPI_ISL_402125.fasta \
    --metadata data/metadata_aligned.tsv


# and plot
python3 scripts/plot_annotated_beast_subtree.py \
    --tree data/weighted_downsampling/19B/19B_location_mcc_median.tre
    --focusTips data/weighted_downsampling/19B/19B_location_mcc_clusters_0.3_0.tsv


# //// DOWNSAMPLE THE GLOBAL SEQUENCES BASED ON MINIMUM DISTANCE ////




# //// 8. DOWNSAMPLE THE GLOBAL SEQUENCES BASED ON SAMPLING DATE ////
python3 scripts/downsample_seqs.py \
    --sequences ${globalSeqs%.fasta}_EHC_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --include ./config/include.tsv \
    --exclude ./config/exclude.tsv \
    --maxDate 2020-03-31 \
    --samplesPerWeek 20 \
    --outName data/temporal_downsampling/ga_focused_aligned_masked_temporal

# subset metadata
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_aligned.tsv \
    > data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv
# get metadata for just GA included sequences
grep "\tGA\t" data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
    > data/temporal_downsampling/ga_focused_aligned_masked_temporal_ga_included_seqs.tsv

# //// 9. BUILD ML TREE OF TEMPORAL DOWNSAMPLING SEQUENCES ////
iqtree2 -redo --polytomy -m TEST -T AUTO -bb 1000 --wbtl -T 4 --prefix data/temporal_downsampling/ga_focused_aligned_masked_temporal -s data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta


# //// 9. ESTIMATE NUMBER OF IMPORTATIONS BASED ON TEMPORAL ML TREE ////
# //// 10. REPEAT ON 100 BOOTSTRAP REPLICATES OF TEMPORAL ML TREE ////
# //// 11. ISOLATE 19B SEQUENCES FROM WEIGHTED ML TREE ////
grep '19B' data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    > data/weighted_downsampling_19B/ga_focused_aligned_masked_weighted_19B_included_seqs.tsv


# //// 12. ADD EHC 19B SEQUENCES TO L84S SEQUENCES AND ALIGNS ////
# todo update this to pull ehc 19B sequences from ML tree
# todo adjust such that we remove EHC seqs from alignment, then concatenate 19B ehc seqs
# todo 1) pull EHC seqs from cluster of beast tree
# 2) remove any of these EHC seqs from the gisaid seqs
# 3) cat cluster seqs and 19B seqs

# for test runs, we didn't have any EHC seqs in teh beast analysis
python3 scripts/align_seqs.py \
    --sequences ${l84sSeqs} \
    --referenceName $refName \
    --maskHead 100 \
    --maskTail 100 \
    --maskSites 11083 15324 21575  \
    --alignType reference \
    --minLength 29000


awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(awk -F'|' '{print $2}' data/weighted_downsampling_old/19B/19B_location_mcc_clusters_0.3_0_match_focal_snps.tsv) \
    data/metadata_aligned.tsv \
    > test.tsv

