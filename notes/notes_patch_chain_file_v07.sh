
# It appears that the start-end coordinates in the target region in the "chain" lines,
# for inverted regions (col 10 = "-"), need to be in the order "lower" "higher".
# Patch the current version:

  cat Wm82v1-Wm82v2_072418.chain | 
    awk -v OFS="\t" '$10~/-/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $12, $11, $13} $10!~/-/ {print}' | 
    awk 'NR==1 {print} $1~/chain/ && NR>1 {print "\n" $0} $1!~/chain/ {print}' \
    > Wm82v1-Wm82v2_072418.patched.chain


# Try CrossMap on the patched file
  CrossMap.py bed Wm82v1-Wm82v2_072418.patched.chain soysnp50k.a1.positions 

  # Error:
			@ 2018-08-05 10:17:25: Read chain_file:  Wm82v1-Wm82v2_072418.patched.chain
			Traceback (most recent call last):
				File "/Users/stevencannon/miniconda3/bin/CrossMap.py", line 1611, in <module>
					(mapTree,targetChromSizes, sourceChromSizes) = read_chain_file(chain_file, print_table = False)
				File "/Users/stevencannon/miniconda3/bin/CrossMap.py", line 274, in read_chain_file
					raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)
			NameError: global name 'header' is not defined

  # It seems the last chain isn't being read correctly:
    chain	9	Gm15	50939160	+	50912071	50913098	scaffold_99	72843	+	69503	70535	23921
    25	1	0
    7	1	0
    25	1	0
    84	1	0
    60	1	0
    826

# Do more stringent filtering on the nucmer results.
# This came from /scratch/abrown1/mummer_unmasked/
# on lathyrus:/scratch/scannon/liftover/
  DIRIN="nucmer_r.v1.q.v2.71218"
  DIROUT="nucmer_r.v1.q.v2.2018-08-05"
  for path in $DIRIN/*.delta; do 
    file=`basename $path .delta`
    delta-filter -q $path -i 99.9 -l 20000 > $DIROUT/$file.filter
    show-coords -rclT $DIROUT/$file.filter > $DIROUT/$file.coords
  done &

# Do mummerplot
  DIROUT="nucmer_r.v1.q.v2.2018-08-05"
  for chr in 01 02 03 04 05 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm01_Gm05/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.2018-08-05"
  for chr in 06 07 08 09 10 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm06_Gm10/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.2018-08-05"
  for chr in 11 12 13 14 15 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm11_Gm15/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.2018-08-05"
  for chr in 16 17 18 19 20 XX; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm16_Gm20/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.2018-08-05"
  for chr in XX; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/scaff/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter
  done &


# Calculate a rough n50 etc. for alignment size for a chromosome (To do: rework into a little script)
  cat nucmer_r.v1.q.v2.2018-08-05/r.v1.q.v2.Gm01.coords | 
    awk 'NR>4 {print $5; sum+=$5} END{print sum}' | 
    sort -nr | awk -v OFS="\t" 'NR==1 {tot=$1} NR>1 {ct++; sum+=$1; print ct, $1, sum/tot}'

# Run modified script on the filtered alignments (length >= 20000 and id >= 99.9(
# from ~/Desktop/liftover/
  ./scripts/delta_to_chain_v05.pl nucmer_r.v1.q.v2.2018-08-05/*filter > nucmer_r.v1.q.v2.2018-08-05.chain

  CrossMap.py bed nucmer_r.v1.q.v2.2018-08-05.chain soysnp50k.a1.positions \
    > nucmer_r.v1.q.v2.2018-08-05.mapped

  # All Gm positions (no scaffolds):
    awk '$1~/^Gm/' nucmer_r.v1.q.v2.2018-08-05.mapped | wc -l
      42795
  # All Gm positions (no scaffolds):
    awk '$1~/^Gm/ && $4!~/Fail/' nucmer_r.v1.q.v2.2018-08-05.mapped | wc -l
      38495
    perl -le 'print 38495/42795'
      0.8995  # too low ... probably.

    Gm01    49225   49225   ->      Gm01    49238   49238
    Gm01    61178   61178   ->      Gm01    61191   61191
    Gm01    81983   81983   Fail
    Gm01    83295   83295   Fail
    Gm01    90419   90419   Fail
    Gm01    92502   92502   Fail
    Gm01    96941   96941   Fail
    Gm01    98352   98352   Fail
    Gm01    138835  138835  Fail
    Gm01    148481  148481  Fail
    Gm01    166873  166873  Fail
    Gm01    183522  183522  Fail
    Gm01    185356  185356  Fail
    Gm01    186315  186315  Fail
    Gm01    225650  225650  Fail
    Gm01    234280  234280  Fail
    Gm01    237375  237375  Fail
    Gm01    243024  243024  Fail
    Gm01    244270  244270  Fail
    Gm01    245175  245175  Fail
    Gm01    281790  281790  Fail
    Gm01    283678  283678  Fail
    Gm01    334575  334575  ->      Gm01    334437  334437
    Gm01    336090  336090  ->      Gm01    335952  335952

  # Check locations of failures relative to the dot plots. Focus on particular chromosomes this way:
    awk '$1~/Gm03/ && ($1==$5 || NF==4)' nucmer_r.v1.q.v2.2018-08-05.mapped
  # or find chromosomes with particularly high failure rates this way:
    awk '$1!~/SCAFF/ && ($1==$5 || NF==4) {print $1 "\t" $5}' nucmer_r.v1.q.v2.2018-08-05.mapped | sort | uniq -c
			 163 Gm01	
			1576 Gm01	Gm01
			 186 Gm02	
			2265 Gm02	Gm02
			 333 Gm03	       high
			1368 Gm03	Gm03
			 148 Gm04	
			1789 Gm04	Gm04
			 132 Gm05	
			1919 Gm05	Gm05
			 168 Gm06	
			1804 Gm06	Gm06
			 299 Gm07	
			1905 Gm07	Gm07
			 182 Gm08	
			2444 Gm08	Gm08
			 186 Gm09	
			1673 Gm09	Gm09
			 198 Gm10	
			2027 Gm10	Gm10
				96 Gm11	      low
			1590 Gm11	Gm11
			 173 Gm12	
			1598 Gm12	Gm12
			 184 Gm13	
			2498 Gm13	Gm13
			 534 Gm14	      high
			1451 Gm14	Gm14
			 269 Gm15	
			2182 Gm15	Gm15
			 166 Gm16	
			1577 Gm16	Gm16
			 148 Gm17	
			1880 Gm17	Gm17
			 333 Gm18	
			2939 Gm18	Gm18
			 244 Gm19	
			2102 Gm19	Gm19
			 158 Gm20	
			1514 Gm20	Gm20

  # The failures appear to correspond with regions of poor alignment, and not to inversions.
  # That's good. The script is likely working.

# Try Annie's prior (more lenient) filtering
  ./scripts/delta_to_chain_v05.pl nucmer_r.v1.q.v2.71218/*filter > nucmer_r.v1.q.v2.71218.chain
  CrossMap.py bed  nucmer_r.v1.q.v2.71218.chain

  CrossMap.py bed nucmer_r.v1.q.v2.71218.chain soysnp50k.a1.positions \
    > nucmer_r.v1.q.v2.71218.mapped

 wc -l *mapped
   42855 nucmer_r.v1.q.v2.2018-08-05.mapped
   43026 nucmer_r.v1.q.v2.71218.mapped

  # All Gm positions (no scaffolds):
    awk '$1~/^Gm/' nucmer_r.v1.q.v2.71218.mapped | wc -l
      42966
  # All Gm positions (no scaffolds):
    awk '$1~/^Gm/ && $4!~/Fail/' nucmer_r.v1.q.v2.71218.mapped | wc -l
      42795
    perl -le 'print 42795/42966'
      0.996  # suspiciously high?

  # check failure rates:
    awk '$1!~/SCAFF/ && ($1==$5 || NF==4) {print $1 "\t" $5}' nucmer_r.v1.q.v2.71218.mapped | sort | uniq -c
				 3 Gm01	
			1743 Gm01	Gm01
				 9 Gm02	
			2440 Gm02	Gm02
				17 Gm03	
			1674 Gm03	Gm03
				 3 Gm04	
			1932 Gm04	Gm04
				12 Gm05	
			2038 Gm05	Gm05
				18 Gm06	
			1955 Gm06	Gm06
				 8 Gm07	
			2205 Gm07	Gm07
				12 Gm08	
			2629 Gm08	Gm08
				 8 Gm09	
			1854 Gm09	Gm09
				 3 Gm10	
			2228 Gm10	Gm10
				 3 Gm11	
			1689 Gm11	Gm11
				 4 Gm12	
			1776 Gm12	Gm12
				 3 Gm13	
			2681 Gm13	Gm13
				 6 Gm14	
			1983 Gm14	Gm14
				 4 Gm15	
			2451 Gm15	Gm15
				 7 Gm16	
			1738 Gm16	Gm16
				 5 Gm17	
			2017 Gm17	Gm17
				17 Gm18	
			3269 Gm18	Gm18
				22 Gm19	
			2329 Gm19	Gm19
				 7 Gm20	
			1667 Gm20	Gm20

# Check plots for the more lenient filtering
# Do mummerplot
  DIROUT="nucmer_r.v1.q.v2.71218"
  for chr in 01 02 03 04 05 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm01_Gm05/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.71218"
  for chr in 06 07 08 09 10 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm06_Gm10/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.71218"
  for chr in 11 12 13 14 15 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm11_Gm15/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.71218"
  for chr in 16 17 18 19 20; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm16_Gm20/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter 
  done &

  DIROUT="nucmer_r.v1.q.v2.71218"
  for chr in XX; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/scaff/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter
  done &

# Conclusion after looking at dot plots and at nucmer_r.v1.q.v2.71218.mapped : 
# Some of the mappings could be spurious, but most look OK. Some post-mapping filtering could be done
# to flag coordinates that are way outside the neighborhood for surrounding markers.

# To do:
# Add usage information etc. to the script
# Make a set of test data
# Test MUMMer filter parameters

#####
# After filtering examination of coords files: use identity >= 99.5 and length >= 1000

# This came from /scratch/abrown1/mummer_unmasked/
# on lathyrus:/scratch/scannon/liftover/
  DIRIN="nucmer_r.v1.q.v2.71218"
  DIROUT="nucmer_r.v1.q.v2.id995_len1k"
  mkdir $DIROUT
  for path in $DIRIN/*.delta; do 
    file=`basename $path .delta`
    delta-filter -q $path -i 99.5 -l 1000 > $DIROUT/$file.filter
    show-coords -rclT $DIROUT/$file.filter > $DIROUT/$file.coords
  done &

# Do mummerplot
  DIROUT="nucmer_r.v1.q.v2.id995_len1k"
  for chr in 01 02 03 04 05 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm01_Gm05/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter &
  done 

  DIROUT="nucmer_r.v1.q.v2.id995_len1k"
  for chr in 06 07 08 09 10 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm06_Gm10/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter &
  done 

  DIROUT="nucmer_r.v1.q.v2.id995_len1k"
  for chr in 11 12 13 14 15 ; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm11_Gm15/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter &
  done 

  DIROUT="nucmer_r.v1.q.v2.id995_len1k"
  for chr in 16 17 18 19 20 XX; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/Gm16_Gm20/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter &
  done 

  DIROUT="nucmer_r.v1.q.v2.id995_len1k"
  for chr in XX; do
    mummerplot -R /scratch/abrown1/mummer_unmasked/Genome/Gmax_v1.1_189.fa \
               -Q /scratch/abrown1/mummer_unmasked/01_pchr_v2/scaff/Gm$chr \
               --png -medium -prefix $DIROUT/Gm$chr $DIROUT/r.v1.q.v2.Gm$chr.filter &
  done 


# Calculate the chain file
  ./scripts/delta_to_chain_v05.pl nucmer_r.v1.q.v2.id995_len1k/*filter > nucmer_r.v1.q.v2.id995_len1k.chain

# Calculate the mapping with CrossMap
  CrossMap.py bed nucmer_r.v1.q.v2.id995_len1k.chain soysnp50k.a1.positions \
    > nucmer_r.v1.q.v2.id995_len1k.mapped

 wc -l *mapped
   42855 nucmer_r.v1.q.v2.2018-08-05.mapped
   43026 nucmer_r.v1.q.v2.id995_len1k.mapped
   43009 nucmer_r.v1.q.v2.id995_len1k.mapped

  # All Gm positions (no scaffolds):
    awk '$1~/^Gm/' nucmer_r.v1.q.v2.id995_len1k.mapped | wc -l
      42949
  # All mapped (not failing) Gm positions (no scaffolds):
    awk '$1~/^Gm/ && $4!~/Fail/' nucmer_r.v1.q.v2.id995_len1k.mapped | wc -l
      42598
    perl -le 'print 42598/42949'
      0.9918  # Probably good

  # check failure rates:
    awk '$1!~/SCAFF/ && ($1==$5 || NF==4) {print $1 "\t" $5}' nucmer_r.v1.q.v2.id995_len1k.mapped | sort | uniq -c
				11 Gm01	
			1735 Gm01	Gm01
				15 Gm02	
			2435 Gm02	Gm02
				49 Gm03	
			1643 Gm03	Gm03
				13 Gm04	
			1922 Gm04	Gm04
				16 Gm05	
			2034 Gm05	Gm05
				26 Gm06	
			1945 Gm06	Gm06
				21 Gm07	
			2192 Gm07	Gm07
				16 Gm08	
			2625 Gm08	Gm08
				15 Gm09	
			1847 Gm09	Gm09
				 9 Gm10	
			2222 Gm10	Gm10
				 5 Gm11	
			1689 Gm11	Gm11
				12 Gm12	
			1768 Gm12	Gm12
				 7 Gm13	
			2676 Gm13	Gm13
				32 Gm14	
			1957 Gm14	Gm14
				11 Gm15	
			2445 Gm15	Gm15
				14 Gm16	
			1730 Gm16	Gm16
				 8 Gm17	
			2013 Gm17	Gm17
				32 Gm18	
			3254 Gm18	Gm18
				29 Gm19	
			2318 Gm19	Gm19
				10 Gm20	
			1663 Gm20	Gm20


# For QC, add temporary names to the SNPs
  cat soysnp50k.a1.positions | awk -v OFS="\t" '{print $1, $2, $3, "Wm82.a1." $1 "." $2}' \
    > soysnp50k.a1.pos_named

# Calculate the chain file
  ./scripts/delta_to_chain.pl nucmer_r.v1.q.v2.id995_len1k/*filter > nucmer_r.v1.q.v2.id995_len1k.chain
  wc nucmer_r.v1.q.v2.id995_len1k.chain

# Calculate the mapping with CrossMap
  CrossMap.py bed sample_out/nucmer_r.v1.q.v2.id995_len1k.chain sample_snps/SNP50K_Wm82.a1.bed \
    > sample_out/SNP50K_Wm82.a2.id995_len1k.bed


#####
# NOTE: See README with workflow at ~/proj/18/crossmap-workflow
# and at https://github.com/soybase/crossmap-workflow

#####
# Evaluate mappings for the SNP50k data for Wm81.a1 and Wm82.a2.
# I received a spreadsheet from David Grant on 2018-08-09, with positions dumped from SoyBase:
#   SoySNP50K_Wm82.a1_vs_Wm82.a2.xls
# I believe the a2 mappings were based on alignments, with positions determined by Nathan
# using a ksh script.

# Join the results for .a1 and .a2 
  cd ~/proj/18/crossmap-workflow/sample_snps

	sort -k4,4 SNP50K_Wm82.a1.bed > SNP50K_Wm82.a1.bed.s
	sort -k4,4 SNP50K_Wm82.a2.bed > SNP50K_Wm82.a2.bed.s

  join -a1 -1 4 -2 4 SNP50K_Wm82.a1.bed.s SNP50K_Wm82.a2.bed.s | 
    awk -v OFS="\t" 'NF==7 {print $1, $2, $3, $4, $5, $6, $7} NF==4 {print $4, $1, $2, $3, "\t\t\t"}' \
    > SNP50K_Wm82.a1.to.Wm82.a2.tsv


  cd ~/proj/18/crossmap-workflow/

  # Note: for the 98 markers with split mappings, pick the first one. This may not always be "correct."
  cat sample_out/SNP50K_Wm82.a2.id995_len1k.bed | 
    awk -v OFS="\t" 'NF==9 {print $6, $7, $8, $9}' | sort -k4,4 | uniq |
      awk '$4==prev && ct<1 {print; ct++} $4!=prev {print; ct=1; prev=$4}' \
    > crossmapped_to_a2.bed
  
  cat sample_snps/SNP50K_Wm82.a2.bed.s | sort -k4,4 | uniq > seq_placed_to_a2.bed

  join -a1 -1 4 -2 4 crossmapped_to_a2.bed seq_placed_to_a2.bed | 
    awk -v OFS="\t" 'NF==7 {print $1, $2, $3, $4, $5, $6, $7} NF==4 {print $4, $1, $2, $3}' \
    > compare_crossmapped_a2.x.seq_placed_a2.tsv

  cat compare_crossmapped_a2.x.seq_placed_a2.tsv |
    awk 'NF==7 {mapped++} 
         NF==4 {unmapped++} 
         END{print "mapped: " mapped; print "unmapped: " unmapped; print "pct mapped: " mapped/(unmapped+mapped)}'
      # mapped: 59648
      # unmapped: 414
      # pct mapped: 0.993107

  cat compare_crossmapped_a2.x.seq_placed_a2.tsv |
    awk '$3==$6 {same_nt++} 
         $3!=$6 {diff_nt++} 
         $2==$5 {same_chr++} 
         $2!=$5 {diff_chr++} 
         END{print "pct same nt:  " same_nt/(same_nt+diff_nt); 
             print "pct same chr: " same_chr/(same_chr+diff_chr)}'
      # pct same nt:  0.6093
      # pct same chr: 0.9930

# Check markers with "close" nucleotide position:
  cat compare_crossmapped_a2.x.seq_placed_a2.tsv |
    awk 'sqrt(($3-$6)^2)<5000 {close_nt++} 
         sqrt(($3-$6)^2)>=5000 {diff_nt++} 
         $2==$5 {same_chr++} 
         $2!=$5 {diff_chr++} 
         END{print "pct close nt:  " close_nt/(close_nt+diff_nt); 
             print "pct same chr: " same_chr/(same_chr+diff_chr)}'
      # pct close nt: 0.9137
      # pct same chr: 0.9930

# Check chromosomes of markers without identical position:
  cat compare_crossmapped_a2.x.seq_placed_a2.tsv |
    awk -v OFS="\t" '$2==$5 && $3!=$6 && sqrt(($3-$6)^2)<5000 {print "DIFF", $2} 
                     $3==$6 {print "SAME", $2}' |
    sort -k2,2 -k1,1 | uniq -c | grep DIFF | awk '$1>1 {print $3 "\t" $1}'
      Gm01	1702
      Gm03	1252
      Gm04	878
      Gm06	2162
      Gm07	1557
      Gm08	3708
      Gm09	2451
      Gm17	139
      Gm19	2305
      Gm20	2122

Spot checking a few …
In all but the first case below, the Crossmap coordinates are the same as the coordinates in GBrowse a2, and the spreadsheet is different – which makes me think that the spreadsheet has NCBI coordinates? – Or at least coordinates different from those in GBrowse.

ss715580806
Crossmap a2:         Gm01   7757883       A
Spreadsheet a2:      Gm01   7758327       B
SoyBase Gbrowse:     Gm01   7758327       B

ss715580734
Crossmap a2:         Gm01   56830220      A      
Spreadsheet a2:      Gm01   56829991      B
SoyBase Gbrowse:     Gm01   56830220      A

ss715602531
Crossmap a2:         Gm01   47796376      A      
Spreadsheet a2:      Gm01   47795458      B
SoyBase Gbrowse:     Gm01   47796376      A

ss715636898
Crossmap a2:         Gm20   21388959      A      
Spreadsheet a2:      Gm20   21388358      B
SoyBase Gbrowse:     Gm20   21388959      A

ss715635164
Crossmap a2:         Gm20   42207155      A      
Spreadsheet a2:      Gm20   42206743      B
SoyBase Gbrowse:     Gm20   42207155      A

ss715602030
Crossmap a2:         Gm08   42884312      A      
Spreadsheet a2:      Gm08   42883394      B
SoyBase Gbrowse:     Gm08   42884312      A
Searching for the 121-base flanking region indicated in the spreadsheet: Gm08:42,883,334..42,883,454 … the SNP is not contained (it is 918 bases away from the center).


#####
# Evaluate mappings for the SNP50k data for Wm81.a1 and Wm82.a2 - WITH A2 DERIVED FROM 
#   the soysnp50k.gff3 from the SoyBase browser.
# I downloaded soysnp50k.gff3 from the a2 browser on 2018-08-10 and renamed it:
  cd ~/proj/18/crossmap-workflow
  mv soysnp50k.gff3.txt soysnp50k_Wm82.a2.gff3

# Convert GFF to minimal BED, sorted by ID
  cat soysnp50k_Wm82.a2.gff3 | awk -v OFS="\t" '$1!~/^##/ {print $1, $4, $5, $9}' |
    perl -pe 's/Name=([^;]+);.+/$1/' | sort -k4,4 > soysnp50k_Wm82.a2.bed



  join -a1 -1 4 -2 4 crossmapped_to_a2.bed soysnp50k_Wm82.a2.bed | 
    awk -v OFS="\t" 'NF==7 {print $1, $2, $3, $4, $5, $6, $7} NF==4 {print $4, $1, $2, $3}' \
    > compare_crossmapped_a2.x.gbrowse_a2.tsv

  cat compare_crossmapped_a2.x.gbrowse_a2.tsv |
    awk 'NF==7 {mapped++} 
         NF==4 {unmapped++} 
         END{print "mapped: " mapped; print "unmapped: " unmapped; print "pct mapped: " mapped/(unmapped+mapped)}'
      # mapped: 59912
      # unmapped: 150
      # pct mapped: 0.997503

  cat compare_crossmapped_a2.x.gbrowse_a2.tsv |
    awk '$3==$6 {same_nt++} 
         $3!=$6 {diff_nt++} 
         $2==$5 {same_chr++} 
         $2!=$5 {diff_chr++} 
         END{print "pct same nt:  " same_nt/(same_nt+diff_nt); 
             print "pct same chr: " same_chr/(same_chr+diff_chr)}'
      # pct same nt:  0.9141
      # pct same chr: 0.9974


# Check markers with "close" nucleotide position:
  cat compare_crossmapped_a2.x.gbrowse_a2.tsv |
    awk 'sqrt(($3-$6)^2)<5000 {close_nt++} 
         sqrt(($3-$6)^2)>=5000 {diff_nt++} 
         $2==$5 {same_chr++} 
         $2!=$5 {diff_chr++} 
         END{print "pct close nt: " close_nt/(close_nt+diff_nt); 
             print "pct same chr: " same_chr/(same_chr+diff_chr)}'
      # pct close nt: 0.9144
      # pct same chr: 0.9974

# Check chromosomes of markers without identical position:
  cat compare_crossmapped_a2.x.gbrowse_a2.tsv |
    awk -v OFS="\t" '$2==$5 && $3!=$6 && sqrt(($3-$6)^2)<5000 {print "CLOSE", $2} 
                     $2==$5 && $3!=$6 && sqrt(($3-$6)^2)>=5000 {print "FAR", $2} 
                     $3==$6 {print "SAME", $2}' |
    sort -k2,2 -k1,1 | uniq -c | grep -v SAME

         1 CLOSE	Gm01
        29 FAR	Gm01
         1 CLOSE	Gm02
        99 FAR	Gm02
       128 FAR	Gm03
       154 FAR	Gm04
         1 CLOSE	Gm05
       975 FAR	Gm05
        65 FAR	Gm06
       517 FAR	Gm07
         1 CLOSE	Gm08
        63 FAR	Gm08
        42 FAR	Gm09
       214 FAR	Gm10
       379 FAR	Gm11
        27 FAR	Gm12
         1 CLOSE	Gm13
      1172 FAR	Gm13
         1 CLOSE	Gm14
       108 FAR	Gm14
        39 FAR	Gm15
       138 FAR	Gm16
         9 FAR	Gm17
       173 FAR	Gm18
         4 CLOSE	Gm19
       240 FAR	Gm19
       189 FAR	Gm20
         1 FAR	scaffold_101
         2 FAR	scaffold_105
       156 FAR	scaffold_21
         6 FAR	scaffold_22
        10 FAR	scaffold_24
         1 CLOSE	scaffold_27
        40 FAR	scaffold_28
         1 FAR	scaffold_31
         1 FAR	scaffold_36
         1 CLOSE	scaffold_41
         2 FAR	scaffold_41
         1 CLOSE	scaffold_72
         1 FAR	scaffold_72
         1 FAR	scaffold_73
         1 CLOSE	scaffold_74
         2 FAR	scaffold_75
         1 FAR	scaffold_91
         2 FAR	scaffold_93
         1 FAR	scaffold_95

# To do: figure out why so many (9%) markers are distant when comparing crossmapping and GBrowse results
  cat compare_crossmapped_a2.x.gbrowse_a2.tsv | awk '$3!=$6 {print $3-$6 "\t" $0}'

# The problem was that the script needed a different transformation for inverted regions:
# transform by 2*(THIS_VALUE-inv_start)-inv_end

### Calculate the chain file:
  scripts/delta_to_chain.pl sample_nucmer/*filter \
    > sample_out/nucmer_r.v1.q.v2.id995_len1k.chain

### Calculate the mapping with CrossMap

  CrossMap.py bed sample_out/nucmer_r.v1.q.v2.id995_len1k.chain sample_snps/SNP50K_Wm82.a1.bed \
    > sample_out/SNP50K_Wm82.a2.id995_len1k.bed

##################################################
### Test with Gm18
  grep Gm18 sample_snps/SNP50K_Wm82.a1.bed > test.a1.Gm18.bed
  
# Check characteristics of the nucmer COORDS data. Note inversions (col 2 below):
  cat sample_nucmer/r.v1.q.v2.Gm18.coords | awk -v OFS="\t" '$5>100000 {print $2-$1, $4-$3, $1, $2, $3, $4 }'

# Check characteristics of the nucmer FILTER data. Note inversions (col 2 below):
  cat sample_nucmer/r.v1.q.v2.Gm18.filter | awk -v OFS="\t" 'NF==7 && $2-$1>100000 {print $2-$1, $4-$3, $0}'
  # An inversion is indicated by query_start > query_end: $3 > $4
  # First large inversion starts at line:
  #   27824914 28003820 31195800 31016909 170 170 0  
  # ... then line 5331 :
  #   28001799 28124274 31018616 30896149 78 78 0
28001799 28124274 31018616 30896149 78 78 0
37
2
9
1
1
1
1
22
97
208
-101
-215
-84
121621
-7
3
0

Chain:
chain 78  Gm18  62308140  + 28001799  28124274  Gm18  58018742  - 30896149  31018616  445
36  1 0
1 1 0
8 5 0
21  1 0
96  1 0
207 1 0
100 0 1
214 0 1
83  0 1
121620  1 0
6 0 1
2 1 0
69

# Nucmer format: http://mummer.sourceforge.net/manual/#nucmeroutput
  A = ABCDACBDCAC$
  B = BCCDACDCAC$
  Delta = (1, -3, 4, 0)
  A = ABC.DACBDCAC$
  B = .BCCDAC.DCAC$

# Get GBrowse file for comparison/validation
  # awk '$1~/^Gm18/' seq_placed_to_a2.bed | sort -k4,4 > test.Gm18.gbrowse.bed

    Gm18  10  150 test1
    Gm18  12  148 test2
    Gm18  13  147 test3
    Gm18  52  110 test4
    Gm18  126 126 test6

##########
# Generate chain and then CrossMap mapping
  #scripts/delta_to_chain.pl sample_nucmer/r.v1.q.v2.Gm18.filter > test.Gm18.chain
  #scripts/delta_to_chain.pl sample_nucmer/r.v1.q.v2.Gm18.filter | perl -pe 's/\tXX.+//' > test.Gm18.chain
  scripts/delta_to_chain.pl test.Gm18.filter | perl -pe 's/\tXX.+//' > test.Gm18.chain

  scripts/delta_to_chain.pl test.Gm18.filter | perl -pe 's/\tXX.+//' |
    awk 'NF==13 {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}
         NF==3 {print $1 "\t" $2 "\t" $3}
         NF==1 {print $1}
         NF==0 {print}' > test.Gm18.chain

  #CrossMap.py bed test.Gm01.chain test.a1.Gm18.bed > test.a2.CM.Gm18.bed

  # In the following command, I have patched CrossMap.py to add "+" or "-" to the output
  CrossMap.py bed test.Gm18.chain test.a1.Gm18.bed | grep -v Fail | grep -v split |     
    awk -v OFS="\t" '{print $1, $2, $7, $10, $4}' | sort -k5,5  > test.a1_to_a2.Gm18.CM_out


  join -a1 -1 5 -2 4 test.a1_to_a2.Gm18.CM_out test.Gm18.gbrowse.bed | 
    awk -v OFS="\t" 'NF==8 {print $7-$4, $1, $2, $3, $4, $5, $6, $7} NF==4 {print "XX", $4, $1, $2, $3}' \
    > compare_crossmapped_Gm18_a2.x.seq_placed_a2.tsv

  # Most map identially to GBrowse v2:
    awk '$1==0' compare_crossmapped_Gm18_a2.x.seq_placed_a2.tsv | wc -l
      4231
    awk '$1!=0' compare_crossmapped_Gm18_a2.x.seq_placed_a2.tsv | wc -l
       177
    # but not those with negative orientation:
    awk '$1!=0 && $6=="-"' compare_crossmapped_Gm18_a2.x.seq_placed_a2.tsv | wc -l
     172
      # These are mostly in the range of -10 million to +10 million; max abs = -58018722

#####
# After script changes to generate FWD and REV chain output:
  scripts/delta_to_chain.pl -fwd test.Gm18.FWD.chain -rev test.Gm18.REV.chain test.Gm18.filter

  # Forward mapping
  CrossMap.py bed test.Gm18.FWD.chain test.a1.Gm18.bed | grep -v Fail | grep -v split |     
    awk -v OFS="\t" '{print $1, $2, $7, $4}' | sort -k4,4  > test.a1_to_a2.Gm18.CM_out_FWD

  join -a1 -1 4 -2 4 test.a1_to_a2.Gm18.CM_out_FWD test.Gm18.gbrowse.bed | 
    awk -v OFS="\t" 'NF==7 {print $7-$4, $1, $2, $3, $4, $5, $6, $7} NF==4 {print "XX", $4, $1, $2, $3}' \
    > compare_crossmapped_Gm18_a2.x.seq_placed_a2_FWD.tsv

  # Reverse mapping: subtract Q positions from the chromosome size.
  chr_size=`head -1 test.Gm18.REV.chain | cut -f9`
  CrossMap.py bed test.Gm18.REV.chain test.a1.Gm18.bed | grep -v Fail | grep -v split |     
    awk -v OFS="\t" -v CHR=$chr_size '{print $1, $2, CHR-$7, $4}' | 
    sort -k4,4 > test.a1_to_a2.Gm18.CM_out_REV

  join -a1 -1 4 -2 4 test.a1_to_a2.Gm18.CM_out_REV test.Gm18.gbrowse.bed | 
    awk -v OFS="\t" 'NF==7 {print $7-$4, $1, $2, $3, $4, $5, $6, $7} NF==4 {print "XX", $4, $1, $2, $3}' \
    > compare_crossmapped_Gm18_a2.x.seq_placed_a2_REV.tsv

  #IT WORKS!  (2018-08-16)

#####
# Run it again clean, with script to recover REV mapping
  scripts/crossmap_delta_to_chain.pl -fwd test.Gm18.FWD.chain -rev test.Gm18.REV.chain test.Gm18.filter

  # Make file of query chromosome & scaffold sizes (derived from filter files)
  grep -h '>' sample_nucmer/*filter | sed 's/>//' | awk -v OFS="\t" '{print $2, $4}' | 
    sort -u > Q_chr_sizes.tsv

  # Forward mapping
  CrossMap.py bed test.Gm18.FWD.chain test.a1.Gm18.bed | grep -v Fail > test.a1.Gm18.map_FWD_TMP
  CrossMap.py bed test.Gm18.REV.chain test.a1.Gm18.bed | grep -v Fail > test.a1.Gm18.map_REV_TMP

  # Recover (translate) the reverse mappings
  ./scripts/crossmap_recover_rev_coords.pl -q Q_chr_sizes.tsv -m test.a1.Gm18.map_REV_TMP



I AM HERE1

  # Try liftover
    liftOver test.a1.Gm18.bed test.Gm18.chain test.a1_to_a2.Gm18.LO_out test.a1_to_a2.Gm18.LO_out_unmap
    # This runs, but no features get mapped.

  # try gff format 
    cat test.a1.Gm18.bed | 
      awk -v OFS="\t" '{print $1, "nucmer", "SNP", $2, $3, 100, "+", ".", "Name=" $4}' > test.a1.Gm18.gff

    liftOver -gff test.a1.Gm18.gff test.Gm18.chain test.a1_to_a2.Gm18.LO_out test.a1_to_a2.Gm18.LO_out_unmap
    # This runs OK.

    cat test.a1_to_a2.Gm18.LO_out | awk -v OFS="\t" '{print $1, $4, $5, $9}' | 
      perl -pe 's/Name=//' | sort -k4,4 > test.a1_to_a2.Gm18.LO_out.bed

    # There are the same number of lines in the liftOver and the crossMap results, and the mappings
    # are identical for the forward alignments and off by one in the reverse alignments:
      paste test.a1_to_a2.Gm18.CM_out test.a1_to_a2.Gm18.LO_out.bed | awk '{print $7-$3}' | sort | uniq -c
        4239 0
         173 1

  # Check alignment 29933547 - 30215942 in the filter(delta) file to see how SNPs are mapping
      29933547 30215942 29039769 28757371 31 31 0
      -245205
      -1
      -1
      0
    cat compare_crossmapped_Gm18_a2.x.seq_placed_a2.tsv | awk '$4>=29933547 && $4<=30215942'
      -221602	ss715630179	Gm18	29954068	29240850	-	Gm18	29019248
      -221602	ss715630180	Gm18	29989702	29205216	-	Gm18	28983614
      -221602	ss715630181	Gm18	30038441	29156477	-	Gm18	28934875
      -221602	ss715630182	Gm18	30063657	29131261	-	Gm18	28909659
      -221602	ss715630183	Gm18	30101339	29093579	-	Gm18	28871977
      -221602	ss715630184	Gm18	30135319	29059599	-	Gm18	28837997
      -221602	ss715630185	Gm18	30172817	29022101	-	Gm18	28800499
      -221602	ss715630186	Gm18	30205453	28989462	-	Gm18	28767860

    # This is the difference between the GBrowse a2 position and the CrossMapped position:
      perl -le 'print 29019248-29240850'
        -221602

      perl -le 'print 58018742-(28757371-29039769)'

    # Try this transformation: 2*(THIS_VALUE-inv_start)-inv_end
      perl -le 'print 2*(29954068-29933547)-30215942'
        -30174900
      perl -le 'print 30215942-29933547'
         282395
      perl -le 'print 29954068-29240850'
         713218

  # Check alignment 31707640 31947636 27273308 27033277 431 431 0

    cat compare_crossmapped_Gm18_a2.x.seq_placed_a2.tsv | awk '$4>=31707640 && $4<=31947636'
      -3712157	ss715630236	Gm18	31712041	30981065	-	Gm18	27268908
      -3712157	ss715630237	Gm18	31736281	30956825	-	Gm18	27244668
      -3712157	ss715630239	Gm18	31777195	30915911	-	Gm18	27203754
      -3712157	ss715630240	Gm18	31805115	30887991	-	Gm18	27175834
      -3712157	ss715630241	Gm18	31838260	30854808	-	Gm18	27142651
      -3712157	ss715630243	Gm18	31882094	30810978	-	Gm18	27098821
      -3712157	ss715630244	Gm18	31915749	30777318	-	Gm18	27065161

    # This is the difference between the GBrowse a2 position and the CrossMapped position:
      perl -le 'print 27268908-30981065'
        -3712157



I AM HERE2
# There seems to be a bug in CrossMap and liftOver in inverted regions.


# Try CrossMap on gff format
  
  cat test.a1.Gm18.bed | 
    awk -v OFS="\t" '{print $1, "nucmer", "SNP", $2, $3, 100, "+", ".", "Name=" $4}' > test.a1.Gm18.gff

  CrossMap.py bed test.Gm18.chain test.a1.Gm18.bed | grep -v Fail | grep -v split |     
    awk -v OFS="\t" '{print $1, $2, $7, $4}' | sort -k4,4  > test.a1_to_a2_gff.Gm18.CM_out
  # Same result as for bed format:
    Gm18	10	58018732	test1
    Gm18	50	58018694	test2
    Gm18	1101	1103	test3
    Gm18	1123	1125	test4
  # Expected result:
    Gm18	10	152	test1
    Gm18	50	110 test2
    Gm18	1101	1103	test3
    Gm18	1123	1125	test4

  # Input chain data (the first two blocks are contrived for testing purposes)
    head -16 test.Gm18.chain
    chain	4	Gm18	62308140	+	1	163	Gm18	58018742	-	1	161	1
    16	2	0
    60	0	4
    10	4	0
    70

    chain	4	Gm18	62308140	+	1001	1161	Gm18	58018742	+	1001	1163	2
    16	0	2
    60	4	0
    10	0	4
    70

    chain	21	Gm18	62308140	+	126611	203483	Gm18	58018742	+	126624	203517	3
    66894	0	21
    9978

#####
# Try liftOver
#   usage:
#      liftOver oldFile map.chain newFile unMapped
#   oldFile and newFile are in bed format by default, but can be in GFF and
#   maybe eventually others with the appropriate flags below.
#   The map.chain file has the old genome as the target and the new genome
#   as the query.
  
    liftOver test.a1.Gm18.bed test.Gm18.chain test.a1_to_a2.Gm18.LO_out test.a1_to_a2.Gm18.LO_out.unmapped
    liftOver test.a1.Gm18.bed test.Gm18.b.chain test.a1_to_a2.Gm18.LO_out test.a1_to_a2.Gm18.LO_out.unmapped
      # complains:  q end mismatch 31459105 vs 31459081 line 2606 of test.Gm18.b.chain

# Try with human genetic data
   cd ~/proj/18/liftover/test.hg19

  # Runs but nothing gets mapped:
    liftOver test.hg19_point.bed hg19ToHg38.over.chain test.hg19ToHg38_point_out.bed test.hg19ToHg38_point_out.unmap

  # Runs OK
    liftOver test.hg19.bed hg19ToHg38.over.chain test.hg19ToHg38_out.bed test.hg19ToHg38_out.unmap

  # modified - works after fixing test alignment
    liftOver test.hg19.bed hg19ToHg38.over_sc.chain test.hg19ToHg38_out_sc.bed test.hg19ToHg38_out_sc.unmap

    liftOver test.hg19.bed3 hg19ToHg38.over_sc.chain test.hg19ToHg38_out_sc.bed test.hg19ToHg38_out_sc.unmap
      grep test  test.hg19ToHg38_out_sc.*
        test.hg19ToHg38_out_sc.bed:chr1	248956364	248956374	test2
        test.hg19ToHg38_out_sc.bed:chr1	248956313	248956314	test3
        test.hg19ToHg38_out_sc.bed:chr1	1010	1011	test4
        test.hg19ToHg38_out_sc.unmap:chr1	10	20	test1


  # Works:
    CrossMap.py bed hg19ToHg38.over.chain test.hg19_point.bed hg38_point.bed

    CrossMap.py bed test.Gm18.b.chain test.a1.Gm18.bed | grep test
      @ 2018-08-14 10:16:18: Read chain_file:  test.Gm18.b.chain
      Gm18	10	10	test1	-=->	Gm18	248956412	248956412	test1  # NOTE: This is chrSize - newPos
      Gm18	50	50	test2	-=->	Gm18	248956374	248956374	test2  # NOTE: This is chrSize - newPos
      Gm18	1101	1101	test3	-=->	Gm18	1103	1103	test3
      Gm18	1123	1123	test4	-=->	Gm18	1125	1125	test4

   

  # Try GFF for SNP data:
    cat test.hg19_point.bed | 
      awk -v OFS="\t" '{print $1, "madeup", "SNP", $2, $3, 0, "+", ".", Name=$4}' > test.hg19_point.gff

    liftOver -gff test.hg19_point.gff hg19ToHg38.over.chain test.hg19ToHg38_gff.out test.hg19ToHg38_gff.out.unmap







