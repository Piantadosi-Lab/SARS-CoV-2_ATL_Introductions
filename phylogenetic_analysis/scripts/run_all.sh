  # global sequence set includes everything marked
# "complete," "high coverage," and "collection date complete"
# from human host 
# on GISAID as of 2021-03-27
# sampled through 2020-03-31
globalSeqs=data/gisaid_hcov-19_2020_03_31_complete_hc_date.fasta
# GISAID global metadata
# this is the Nextmeta file
# maybe no longer be availabe 
# feel free to email me and I will share if you are a GISAID member
# mmart59@emory.edu
# downloaded on 2021-03-27
globalData=data/metadata.tsv
# EHC sequences
ehcSeqs=data/Early_GA_seqs_51.fasta
# EHC metadata
ehcData=data/Early_GA_seqs_metadata_51.tsv
# TODO update this
l84sSeqs=data/gisaid_hcov-19_2020_03_31_l84s_complete_hc_date.fasta
# reference information
refSeq='data/EPI_ISL_402125.fasta'
refName='EPI_ISL_402125'
# US case data
# from Hopkins Github
usCaseData=data/time_series_covid19_confirmed_US.csv

# sets up folder structure
mkdir data
mkdir config
mkdir data/weighted_downsampling
mkdir data/19B_subclade
mkdir data/l84s
mkdir data/temporal_downsampling
mkdir data/travel_analysis

# --------------------------------------------------------------------------------------#
# //// 1. BASIC FORMATTING OF EHC SEQUENCES AND METADATA ////
# --------------------------------------------------------------------------------------#
# formats the ehc_seqs header so they can be split on a "|"
# also removes "amp" and "merged" from file names
sed '/\>/ s//\>\|/g' $ehcSeqs | \
    sed 's/\_amp//g' | \
    sed 's/\_merged//g' | \
    sed 's/_new//g' | \
    sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' > ${ehcSeqs%.fasta}_format.fasta 
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

# the following were uploaded by the CDC and may be mixed up wrt sequence/accession
echo '# uploaded by CDC and may be mixed up with regard to name/sequence
1\tEPI_ISL_2790694
1\tEPI_ISL_2790695
1\tEPI_ISL_2790696
1\tEPI_ISL_2790697
1\tEPI_ISL_2790698
1\tEPI_ISL_2790699
1\tEPI_ISL_2790700
1\tEPI_ISL_2790701
1\tEPI_ISL_2790702
1\tEPI_ISL_2790703
1\tEPI_ISL_2790704
1\tEPI_ISL_2790705
1\tEPI_ISL_2790706
1\tEPI_ISL_2790707
1\tEPI_ISL_2790708
1\tEPI_ISL_2790709
1\tEPI_ISL_2790710
1\tEPI_ISL_2790711
1\tEPI_ISL_2790712
1\tEPI_ISL_2790713
1\tEPI_ISL_2790714
1\tEPI_ISL_2790715' \
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


# --------------------------------------------------------------------------------------#
# //// 2. GENERATE A SINGLE FORMATTED METADATA FILE FOR BOTH EHC AND GISAID DATA
# --------------------------------------------------------------------------------------#
# human samples only
# remove samples with travel history (sample location doesn't match exposure location)
# remove ambiguous dates
# remove cruise related samples
# exclude EHC sequences from GISAID data so it's not double counted
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

# location columns for phylogenetics
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
    --minLength 28000


# subset the metadata to just the aligned and filtered sequences 
awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
    <(grep ">" ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta |\
    sed 's/\>//g' | awk -F'|' '{print $2}') \
    ${globalData%.tsv}_combined.tsv \
    > ${globalData%.tsv}_aligned.tsv


# --------------------------------------------------------------------------------------#
# //// 3. CALCULATE DOWNSAMPLING WEIGHTS ////
# --------------------------------------------------------------------------------------#
# downsampling scheme is two tiered: 
# 1) weightedd on a per country basis related to the number of 
# cases in each country by March 31st
# 2) 1/(D+1) where D is the maximum distance 
# between that sequence and all of the GA sequences
# // 3A. COUNTRY CASE WEIGHTS //
python3 scripts/calc_case_weights.py \
    --globalDat data/time_series_covid19_confirmed_global.csv \
    --usDat data/time_series_covid19_confirmed_US.csv \
    --dateMax 03/31/20 \
    --focalStates Georgia


# // 3B. DISTANCE WEIGHTS //
# first get all GA sequences from global alignment
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

# --------------------------------------------------------------------------------------#
# //// 4. DOWNSAMPLE THE GLOBAL SEQUENCES BASED ON MINIMUM DISTANCE ////
# --------------------------------------------------------------------------------------#
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

# get fasta of just GA included sequences
python3 scripts/get_seqs.py \
    --seqs data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --getSeqs data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv \
    --getSeqsCol 1 \
    --getSeqsField 0 \
    --outName data/weighted_downsampling/ga_focused_aligned_masked_weighted_GA

# todo refresh this
# pangolin web app used to assign pango lineages on August 26, 2021
# add date col to pangolin 

awk -F'\t' 'FNR==NR {date[$2] = $3; next} {split($1, columns, ","); split(columns[1], name, "|"); print name[2]"\t"date[name[2]]"\t"columns[2]}' \
    data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv \
    <(tail -n +2 data/weighted_downsampling/ga_focused_aligned_masked_weighted_GA_pangolin.csv) \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_GA_pangolin_dates.csv




# --------------------------------------------------------------------------------------#
# //// 5. BUILD ML TREE OF WEIGHTED DOWNSAMPLED SEQS ////
# --------------------------------------------------------------------------------------#
iqtree2 -redo --polytomy -m TEST -bb 1000 --wbtl -T 4 \
    --prefix data/weighted_downsampling/ga_focused_aligned_masked_weighted \
    -s data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta


# --------------------------------------------------------------------------------------#
# //// 6. ESTIMATE NUMBER OF IMPORTATIONS BASED ON WEIGHTED ML TREE ////
# --------------------------------------------------------------------------------------#
# run aanlysis
python3 scripts/estimate_importations.py \
    --trees data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile \
    --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    --sequences data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
    --outgroup $refName \
    --biasCorrection 2.5 \
    --clockRate 0.001 \
    --regions GA \
    --metadataLocCol 6 \
    --saveAln 

# plot RTT
python3 scripts/plot_rtt.py \
    --data data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_rtt.csv \
    --outName figures/weighted_downsampling_rtt


# --------------------------------------------------------------------------------------#
# //// 7. REPEAT ON 100 BOOTSTRAP REPLICATES OF WEIGHTED ML TREE ////
# --------------------------------------------------------------------------------------#
refName='EPI_ISL_402125'
python3 scripts/estimate_importations.py \
            --trees data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot \
            --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
            --sequences data/weighted_downsampling/ga_focused_aligned_masked_weighted.fasta \
            --outgroup $refName \
            --biasCorrection 2.5 \
            --clockRate 0.001 \
            --regions GA \
            --metadataLocCol 6 \
            --bootstrapReps 0 100

# --------------------------------------------------------------------------------------#
# //// 8. RESULTS FROM STEPS 7 AND 8 ////
# --------------------------------------------------------------------------------------#
# RERUN ANYTHING WHICH USES GA TRAVEL SEQS FILE
python3 scripts/plot_annotated_tree.py \
    --tree data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick \
    --treeStates data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_node_states.csv \
    --treeBootstrapDir "data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv" \
    --config config/annotated_tree.json \
    --metadata data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    --travelData <(awk -F'\t' '{print $2}' data/ga_travel_seqs.tsv) \
    --highlightMRCA data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0_family.tsv

    

# plot ML divergence tree
python3 scripts/plot_divergence_tree.py \
    --tree data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined.newick\
    --config config/annotated_global_divergence_tree.json \
    --metadata data/metadata_combined.tsv \
    --metadataRegionCol 6

# estimate number of introductions
python3 scripts/estimate_n_importations_ml.py \
    --treeBootstrapDir "data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv" 

# filter importations by known travelers
python3 scripts/estimate_n_importations_ml.py\
    --treeBootstrapDir "data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv"  \
    --filterSeqs data/ga_travel_seqs.tsv

# generate table for ML tree
awk -F'\t' 'FNR==NR{split($1, a, "|"); hash[a[2]]; next} \
    {label="True"; if ($2 in hash) label="False"; print $1"\t"$2"\t"label}' \
    data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_bad_seqs.tsv \
    data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv \
    > tables/weighted_downsampling_included.tsv


# --------------------------------------------------------------------------------------#
# //// 9. GET MUTATIONAL PROFILE SHARED BETWEEN GRANDPARENT OF MRCA OF 19B SUBCLADE ////
# --------------------------------------------------------------------------------------#
# cluster GA seqs in tree
python3 scripts/cluster_seqs.py \
    --tree  data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick  \
    --seqNames data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv 
# subset to just cluster 0
# always the biggest cluster, they're sorted
grep '\t0' data/weighted_downsampling/ga_focused_aligned_masked_weighted_clusters_0.3.tsv | \
    awk '{print $1}' \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv
# get the pangolin lineage of these sequences
# todo this
grep -f data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv \
    data/weighted_downsampling/ga_focused_aligned_masked_weighted_GA_pangolin.csv \
    > data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0_pangolin.tsv

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
# //// 10. GET ALL AVAILABLE SEQUENCES WHICH MATCH THIS MUTATIONAL PROFILE ////
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


# --------------------------------------------------------------------------------------#
# //// 11. BUILD ML TREE OF 19B SUBCLADE ////
# --------------------------------------------------------------------------------------#
# build ML tree
# todo check iqtree version
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
# get remaining metadata
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


# second ML plot, this time labeling non-USA tips sampled after March 1st
# all from 2020 so can just look at month
python3 scripts/plot_divergence_tree.py \
    --tree data/19B_subclade/19B_subclade_family_clockfilter.newick \
    --config config/annotated_divergence_tree2.json \
    --metadata <(grep "2020-03" data/19B_subclade/19B_subclade_family.tsv | awk -F'\t' '$5!="USA"{print $0}') \
    --metadataRegionCol 7 \
    --highlightMRCA <(grep "USA" data/19B_subclade/19B_subclade_family.tsv)



 
# --------------------------------------------------------------------------------------#
# //// 12. GENERATE BEAST XML FILE AND RUN////
# --------------------------------------------------------------------------------------#
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



# note spelling mistake in alnName parameter
python3 scripts/generate_xml.py \
    --xmlTemplate config/beast2_trait_template.xml \
    --outName data/19B_subclade/19B_subclade \
    --seqs data/19B_subclade/19B_subclade_family_clockfilter.fasta \
    --metadata data/19B_subclade/19B_subclade_beast_include.tsv \
    --metadataTraitCol 7 \
    --alnName 19B_subclade \
    --traitName location 

# run beast
# todo add command
# generate MCC tree
/Applications/BEAST\ 2.6.3/bin/treeannotator -heights median -b 10  -lowMem \
    data/19B_subclade/19B_subclade19B_subclade_location_tree_with_trait.trees \
    data/19B_subclade/19B_location_mcc.tre


# --------------------------------------------------------------------------------------#
# //// 14. STATISTICS AND PLOT FOR BEAST ANALYSIS ////
# --------------------------------------------------------------------------------------#
# parse set of beast trees
python3 scripts/parse_beast_trees.py \
    --trees data/19B_subclade/19B_subclade19B_subclade_location_tree_with_trait.trees


# gets statistics
# todo make this more efficient!
python3 scripts/estimate_n_importations_beast.py \
    --trees data/19B_subclade/19B_subclade19B_subclade_location_tree_with_trait.pkl \
    --focalRegion GeorgiaUSA

# generate table
python3 scripts/generate_beast_table.py \
    --log data/19B_subclade/19B_subclade19B_subclade.log

# plot tree
python3 scripts/plot_annotated_beast.py \
    --tree data/19B_subclade/19B_location_mcc.tre \
    --config config/annotated_beast_tree.json \
    --treeFocusTips data/weighted_downsampling/ga_focused_aligned_masked_weighted_cluster_0.tsv \
    --nImportations data/19B_subclade/19B_subclade19B_subclade_location_tree_with_trait_importations_count.tsv

# --------------------------------------------------------------------------------------#
# //// 13. GENERATE RANDOMLY DOWNSAMPLED BEAST RUNS ////
# --------------------------------------------------------------------------------------#
python3 scripts/downsample_metadata.py \
    --metadata data/19B_subclade/19B_subclade_beast_include.tsv \
    --groupCol 7 \
    --nPerGroup 5 \
    --nReps 5 \
    --timeGroup week \
    --dateCol 2



for i in 0 $(seq 4 $END); 
do 
python3 scripts/generate_xml.py \
        --xmlTemplate config/beast2_trait_template.xml \
        --outName data/19B_subclade/19B_subclade_sampled_5_rep_${i} \
        --seqs data/19B_subclade/19B_subclade_family_clockfilter.fasta \
        --metadata data/19B_subclade/19B_subclade_beast_include_sampled_5_rep_${i}.tsv \
        --metadataTraitCol 7 \
        --alnName 19B_subclade_sampled_5_rep_${i} \
        --traitName location 
done

# run beast on cluster

# get mcc tree for each and process
for i in 0 $(seq 4 $END); 
do 
/Applications/BEAST\ 2.6.3/bin/treeannotator -heights median -b 10  -lowMem \
    data/19B_subclade/19B_subclade_sampled_5_rep_${i}_location_tree_with_trait.trees \
    > data/19B_subclade/19B_subclade_sampled_5_rep_${i}_location_tree_with_trait_mcc.tr
done


for i in 0 $(seq 4 $END); 
do 
#python3 scripts/parse_beast_trees.py \
#    --trees data/19B_subclade/19B_subclade_sampled_5_rep_${i}_location_tree_with_trait.trees
# gets statistics
# todo make this more efficient!
python3 scripts/estimate_n_importations_beast.py \
    --trees data/19B_subclade/19B_subclade_sampled_5_rep_${i}_location_tree_with_trait.pkl \
    --focalRegion GeorgiaUSA
done

# combine parsed files
rm -rf "data/19B_subclade/19B_subclade_sampled_5_rep_*_location_tree_with_trait_importations_count.tsv"

for i in 0 $(seq 4 $END); 
do 
awk -F'\t' -v rep=$i '{print rep"\t"$0}' data/19B_subclade/19B_subclade_sampled_5_rep_${i}_location_tree_with_trait_importations_count.tsv \
    >> "data/19B_subclade/19B_subclade_sampled_5_rep_*_location_tree_with_trait_importations_count.tsv"

done



# process all downsampled BEAST runs
python3 scripts/process_downsampled_beast.py \
    --trees "data/19B_subclade/19B_subclade_sampled_5_rep_*_mcc.tre" \
    --metadata <(awk -F'\t' '$5=="USA"{print $0}' data/19B_subclade/19B_subclade_beast_include.tsv)

# plot results from downsamples
python3 scripts/plot_downsampled_beast.py \
    --downsampleData data/19B_subclade/19B_subclade_sampled_5_rep_*_mcc_parsed.tsv \
    --config config/downsampled_beast.json \
    --nIntroductions "data/19B_subclade/19B_subclade_sampled_5_rep_*_location_tree_with_trait_importations_count.tsv"

# generate included sequences table
awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $0"\t"label}' \
    data/19B_subclade/19B_subclade_beast_include_sampled_5_rep_4.tsv \
    <(awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $0"\t"label}' \
        data/19B_subclade/19B_subclade_beast_include_sampled_5_rep_3.tsv \
        <(awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $0"\t"label}' \
            data/19B_subclade/19B_subclade_beast_include_sampled_5_rep_2.tsv \
            <(awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $0"\t"label}' \
                data/19B_subclade/19B_subclade_beast_include_sampled_5_rep_1.tsv \
                <(awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $0"\t"label}' \
                    data/19B_subclade/19B_subclade_beast_include_sampled_5_rep_0.tsv \
                    <(awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $0"\t"label}' \
                        data/19B_subclade/19B_subclade_beast_include.tsv \
                        <(awk -F'\t' 'FNR==NR{hash[$2]; next}{label="False"; if ($2 in hash) label="True"; print $1"\t"$2"\t"label}' \
                            data/19B_subclade/19B_subclade_family_clockfilter.tsv  \
                            data/19B_subclade/19B_subclade_family.tsv)))))) \
    > tables/19B_subclade_family.tsv

# generate results table
for i in 0 $(seq 4 $END); 
do 
    python3 scripts/generate_beast_table.py \
    --log  data/19B_subclade/19B_subclade_sampled_5_rep_${i}.log
done


# --------------------------------------------------------------------------------------#
# //// 15. GET "FUTURE" SEQUENCES WHICH MATCH MUTATIONAL PROFILE of GA 19B subclade ////
# --------------------------------------------------------------------------------------#
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

# subset to just first two column for table
awk -F'\t' '{print $1"\t"$2}' \
    data/l84s/l84s_19B_subclade_match_focal_snps.tsv \
    > tables/l84s_19B_subclade_match_focal_snps.tsv

python3 scripts/compare_seqs.py \
    --seqs data/l84s/l84s_filtered_aligned_ref_filtered_masked_noref.fasta \
    --focalSeqs data/l84s/l84s_19B_subclade_earliest.tsv \
    --refSeq data/EPI_ISL_402125.fasta \
    --focalSeqsField 0


# --------------------------------------------------------------------------------------#
# //// 16. PLOT NUMBER OF MATCHING SEQUENCES PER WEEK ////
# --------------------------------------------------------------------------------------#
python3 scripts/plot_19B_subclade.py \
    --mutationalProfile data/weighted_downsampling/19B_subclade_focal_snps_shared.tsv \
    --config config/19B_subclade.json \
    --matchMutationalProfile data/l84s/l84s_19B_subclade_match_focal_snps.tsv \


# --------------------------------------------------------------------------------------#
# //// 17. DOWNSAMPLE THE GLOBAL SEQUENCES BASED ON SAMPLING DATE ////
# --------------------------------------------------------------------------------------#
python3 scripts/downsample_seqs.py \
    --sequences ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --include ./config/include.tsv \
    --exclude ./config/exclude.tsv \
    --maxDate 2020-03-31 \
    --samplesPerWeek 20 \
    --outName data/temporal_downsampling/ga_focused_aligned_masked_temporal

# # get clade assignments for downsampled sequences
nextclade \
    --input-fasta data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta \
    --output-tsv data/temporal_downsampling/ga_focused_aligned_masked_temporal_clades.tsv


# get metadata for all included sequences and combine with clades
awk -F'\t' 'FNR==NR {split($1, b, "|"); a[b[2]]=$2; next}{print $0"\t"a[$2]}' \
    data/temporal_downsampling/ga_focused_aligned_masked_temporal_clades.tsv \
    <(awk -F'\t' 'FNR==NR {hash[$1]; next} $2 in hash' \
        <(grep ">" data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta |\
        sed 's/\>//g' | awk -F'|' '{print $2}') \
        ${globalData%.tsv}_aligned.tsv) \
    > data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv



# get metadata for just GA included sequences
grep "\tGA\t" data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
    > data/temporal_downsampling/ga_focused_aligned_masked_temporal_ga_included_seqs.tsv


# --------------------------------------------------------------------------------------#
# //// 18. BUILD ML TREE OF TEMPORAL DOWNSAMPLING SEQUENCES ////
# --------------------------------------------------------------------------------------#
iqtree2 -redo --polytomy -m TEST -T AUTO -bb 1000 --wbtl -T 4 --prefix data/temporal_downsampling/ga_focused_aligned_masked_temporal -s data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta


# --------------------------------------------------------------------------------------#
# //// 19. ESTIMATE NUMBER OF IMPORTATIONS BASED ON TEMPORAL ML TREE ////
# --------------------------------------------------------------------------------------#
python3 scripts/estimate_importations.py \
            --trees data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile \
            --metadata data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
            --sequences data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta \
            --outgroup $refName \
            --biasCorrection 2.5 \
            --clockRate 0.001 \
            --regions GA \
            --metadataLocCol 6 \
            --saveAln 

# --------------------------------------------------------------------------------------#
# //// 20. REPEAT ON 100 BOOTSTRAP REPLICATES OF TEMPORAL ML TREE ////
# --------------------------------------------------------------------------------------#
refName='EPI_ISL_402125'
python3 scripts/estimate_importations.py \
            --trees data/temporal_downsampling/ga_focused_aligned_masked_temporal.ufboot \
            --metadata data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
            --sequences data/temporal_downsampling/ga_focused_aligned_masked_temporal.fasta \
            --outgroup $refName \
            --biasCorrection 2.5 \
            --clockRate 0.001 \
            --regions GA \
            --metadataLocCol 6 \
            --bootstrapReps 0 100


# --------------------------------------------------------------------------------------#
# //// 21. PLOT TEMPORAL DOWNSAMPLING TREE ////
# --------------------------------------------------------------------------------------#
python3 scripts/plot_annotated_tree.py \
    --tree data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_refined_time.newick \
    --treeStates data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_refined_node_states.csv \
    --treeBootstrapDir "data/temporal_downsampling/ga_focused_aligned_masked_temporal.ufboot_tres/*/*_importations.csv" \
    --config config/annotated_temporal_tree.json \
    --metadata data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
    --travelData <(awk -F'\t' '{print $2}' data/ga_travel_seqs.tsv) 

# plot RTT
python3 scripts/plot_rtt.py \
    --data data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_rtt.csv \
    --outName figures/temporal_downsampling_rtt

# get # of introductions
python3 scripts/estimate_n_importations_ml.py \
    --treeBootstrapDir "data/temporal_downsampling/ga_focused_aligned_masked_temporal.ufboot_tres/*/*_importations.csv" 

# generate table

awk -F'\t' 'FNR==NR{split($1, a, "|"); hash[a[2]]; next} \
    {label="True"; if ($2 in hash) label="False"; print $1"\t"$2"\t"label}' \
    data/temporal_downsampling/ga_focused_aligned_masked_temporal.treefile_tres/0/0_bad_seqs.tsv \
    data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv \
    > tables/temporal_downsampling_included.tsv





# --------------------------------------------------------------------------------------#
# //// 22. GET LIST OF ALL INCLUDED SEQUENCES ////
# --------------------------------------------------------------------------------------#
# todo this is not correct! used all sequences
cat <(awk -F'\t' '{print $2}' data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv | grep -v "GA-EHC") \
    <(awk -F'\t' '{print $2}' data/temporal_downsampling/ga_focused_aligned_masked_temporal_all_included_seqs.tsv | grep -v "GA-EHC") \
    <(awk -F'\t' '{print $2}' data/l84s/l84s_aligned.tsv | grep -v "GA-EHC") \
    <(awk -F'\t' '{print $2}' config/exclude.tsv | grep -v "#"  | grep -v "GA-EHC" | grep -v "^$") | \
    uniq > \
    data/all_included_seqs.txt


# --------------------------------------------------------------------------------------#
# //// 23. NUMBER OF CASES AND SEQUENCES OVER TIME ////
# --------------------------------------------------------------------------------------#
python3 scripts/plot_case_data.py \
    --caseData data/time_series_covid19_confirmed_US.csv \
    --seqData data/weighted_downsampling/ga_focused_aligned_masked_weighted_ga_included_seqs.tsv

# --------------------------------------------------------------------------------------#
# //// 24.  TRAVEL ANALYSIS ////
# --------------------------------------------------------------------------------------#
# compare to travel regions
python3 scripts/calc_distance.py \
    --seqs ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --config <(awk -F'\t' '{print $2"\t"$3}' data/ga_travel_seqs.tsv) \
    --refSeq $refSeq \
    --outDir data/travel_analysis \
    > data/travel_analysis/log_travel.txt

# compare to GA sequences
python3 scripts/calc_distance.py \
    --seqs ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --config <(awk -F'\t' '{print $2"\tGeorgiaUSA"}' data/ga_travel_seqs.tsv) \
    --refSeq $refSeq \
    --outDir data/travel_analysis \
    > data/travel_analysis/log_ga.txt

# now compare to all regions
python3 scripts/calc_distance.py \
    --seqs ${globalSeqs%.fasta}_EHC_metadata_aligned_ref_filtered_masked.fasta \
    --metadata ${globalData%.tsv}_aligned.tsv \
    --config <(awk -F'\t' '{print $2}' data/ga_travel_seqs.tsv) \
    --refSeq $refSeq \
    --outDir data/travel_analysis \
    > data/travel_analysis/log_all.txt

# now make the figures 
Rscript scripts/plot_traveler_dists.R

# and table
python3 scripts/generate_traveler_table.py --inDir data/travel_analysis --metadata data/ga_travel_seqs.tsv