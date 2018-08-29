
#####
# Convert MUMmer delta-format files (.delta or .filter) to chain format. 
# Save results in two output files: one for forward alignments and one for inversions (reverse alignments).
# This can be run on a single concatenated file of MUMmer results:
  cat sample_nucmer/*filter > work/all_r.v1.q.v2.filter
  scripts/crossmap_delta_to_chain.pl -fwd work/all_FWD.chain -rev work/all_REV.chain work/all_r.v1.q.v2.filter

# Use the chain files to map features in a supplied BED file. Do this for both the FWD and REV chain files.
  CrossMap.py bed work/all_FWD.chain sample_snps/SNP50K_Wm82.a1.bed | grep -v Fail > work/all_map_a1_a2_fwd.tsv
  CrossMap.py bed work/all_REV.chain sample_snps/SNP50K_Wm82.a1.bed | grep -v Fail > work/all_map_a1_a2_rev_TMP

# Make file of query chromosome & scaffold sizes (derived from filter files)
  grep -h '>' sample_nucmer/*filter | sed 's/>//' | awk -v OFS="\t" '{print $2, $4}' |
    sort -u > work/Q_chr_sizes.tsv

# Recover (translate) the reverse mappings
  scripts/crossmap_recover_rev_coords.pl -q work/Q_chr_sizes.tsv \
    -m work/all_map_a1_a2_rev_TMP -o work/all_map_a1_a2_rev.tsv

# Combined results:
  cat work/all_map_a1_a2_fwd.tsv work/all_map_a1_a2_rev.tsv > sample_out/all_map_a1_a2.tsv


#####
# Evaluate mappings for the SNP50k data for Wm81.a1 and Wm82.a2 - WITH A2 DERIVED FROM
#   the soysnp50k.gff3 from the SoyBase browser.
# I downloaded soysnp50k.gff3 from the a2 browser on 2018-08-10 and renamed it:
  mv sample_snps/soysnp50k.gff3.txt sample_snps/soysnp50k_Wm82.a2.gff3

# Convert GFF to minimal BED, sorted by ID
  cat sample_snps/soysnp50k_Wm82.a2.gff3 | awk -v OFS="\t" '$1!~/^##/ {print $1, $4, $5, $9}' |
    perl -pe 's/Name=([^;]+);.+/$1/' | sort -k4,4 > sample_snps/soysnp50k_Wm82.a2.bed

  cut -f6,7,8,9 sample_out/all_map_a1_a2.tsv | sort -k4,4 > sample_out/all_map_a2.bed

  join -a1 -1 4 -2 4 sample_out/all_map_a2.bed sample_snps/soysnp50k_Wm82.a2.bed |
    awk -v OFS="\t" 'NF==7 {print $1, $2, $3, $4, $5, $6, $7} NF==4 {print $4, $1, $2, $3}' \
    > sample_out/compare_crossmapped_a2.x.gbrowse_a2.tsv

  cat sample_out/compare_crossmapped_a2.x.gbrowse_a2.tsv |
    awk 'NF==7 {mapped++}
         NF==4 {unmapped++}
         END{print "mapped: " mapped; print "unmapped: " unmapped; print "pct mapped: " mapped/(unmapped+mapped)}'
      # mapped: 60656
      # unmapped: 152
      # pct mapped: 0.9975

  cat sample_out/compare_crossmapped_a2.x.gbrowse_a2.tsv | 
    awk -v OFS="\t" '$3==$6 {same++} $3!=$6 {not++} END{print same/(same+not)}'
      0.996843  # Result: 99.68% of ALL locations are the same (most of the mismatches don't map at all).)

  cat sample_out/compare_crossmapped_a2.x.gbrowse_a2.tsv | 
    awk -v OFS="\t" '$3==$6 && NF==7 {same++} $3!=$6 && NF==7 {not++} END{print same/(same+not)}'
      0.999341  # Result: 99.93% of the MAPPED locations are the same.


