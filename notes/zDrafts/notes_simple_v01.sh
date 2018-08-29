
#####
# Convert MUMmer delta-format files (.delta or .filter) to chain format. 
# Save results in two output files: one for forward alignments and one for inversions (reverse alignments).
# This can be run on a single concatenated file of MUMmer results:
  cat sample_nucmer/*filter > all_r.v1.q.v2.filter
  scripts/crossmap_delta_to_chain.pl -fwd test.Gm18.FWD.chain -rev test.Gm18.REV.chain all_r.v1.q.v2.filter

  # Forward mapping
  CrossMap.py bed test.Gm18.FWD.chain test.a1.Gm18.bed | grep -v Fail > test.a1.Gm18.map_FWD_TMP
  CrossMap.py bed test.Gm18.REV.chain test.a1.Gm18.bed | grep -v Fail > test.a1.Gm18.map_REV_TMP

  # Make file of query chromosome & scaffold sizes (derived from filter files)
  grep -h '>' sample_nucmer/*filter | sed 's/>//' | awk -v OFS="\t" '{print $2, $4}' |
    sort -u > Q_chr_sizes.tsv

  # Recover (translate) the reverse mappings
  ./scripts/crossmap_recover_rev_coords.pl -q Q_chr_sizes.tsv -m test.a1.Gm18.map_REV_TMP

#####
# Run on all nucmer (filter) files
  
  cat sample_nucmer/*filter > all_r.v1.q.v2.filter

  # Generate chain files from MUMmer files
  scripts/crossmap_delta_to_chain.pl -fwd all_FWD.chain -rev all_REV.chain all_r.v1.q.v2.filter

  # Forward mapping
  CrossMap.py bed all_FWD.chain sample_snps/SNP50K_Wm82.a1.bed | grep -v Fail > all_map_FWD.tsv
  CrossMap.py bed all_REV.chain sample_snps/SNP50K_Wm82.a1.bed | grep -v Fail > all_map_REV_TMP

  # Make file of query chromosome & scaffold sizes (derived from filter files)
  grep '>' all_r.v1.q.v2.filter | sed 's/>//' | awk -v OFS="\t" '{print $2, $4}' | sort -u > Q_chr_sizes.tsv
    # NOTE: This misses a few scaffolds (XX but why?). Go to the full assembly to get all:
    zcat glyma.Wm82.gnm2.DTC4.genome_main.fna.gz | seqlen.awk | perl -pe 's/glyma.Wm82.gnm2.//' > Q_chr_sizes.tsv

  # Recover (translate) the reverse mappings
  ./scripts/crossmap_recover_rev_coords.pl -q Q_chr_sizes.tsv -m all_map_REV_TMP -o all_map_REV.tsv


