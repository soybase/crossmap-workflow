# crossmap-workflow
The objective of this workflow is to develop a mapping between genomic assemblies, starting from nucmer comparisons of the assemblies. 
This is facilitated by the script "crossmap_delta_to_chain.pl", which generates a chain file from MUMmer's delta-format files. 
Then, given the chain file, CrossMap (or equivalently, liftOver) can generate a new mapping of genomic features.

## Basic workflow:
Run nucmer on each chromosome and scaffold against the whole genome assembly (to allow each region
to find its best home on the whole genome), then filter using the MUMmer utility delta-filter
(recommend parameters for a lenient mapping between genome assemblies: -i 99.5 -l 1000).

Then use the script crossmap_delta_to_chain.pl in this repository to generate "chain" files from the
delta-format .filter files. 

This process requires a couple of additional steps due to an apparent problem with liftOver and CrossMap. As of late 2018, it seems that liftOver and CrossMap are doing the wrong thing
with inversions, so the crossmap_delta_to_chain.pl script produces two chain files: one for forward
alignments and one for inversions. The chain file for inversions transforms the
coords as the query chromosome length minus the original coordinates, in + orientation.

Then use the CrossMap.py script (http://crossmap.sourceforge.net) to generate a mapping of a file of genomic coordinates. 

Finally, recover the corrected coordinates from the mapped features in inverted regions, using the script crossmap_recover_rev_coords.pl.

The instructions below assume local installation of the MUMmer package and CrossMap. Note that CrossMap can be easily installed using Conda, e.g.
```
  conda install -c bioconda mummer 
  conda install -c bioconda crossmap
```

## Example using soybean data
The following example maps a set of SNP positions, in simplified BED format, between the soybean (Glycine max) Williams 82 assembly 1 (2010 JGI assembly version) and Williams 82 assembly 2.

### Filter the raw nucmer output:
Note: see the RESULTS in sample_nucmer/. To save space, the raw delta files have not been included.
```
  for path in sample_nucmer/*.delta; do
    file=`basename $path .delta`
    delta-filter -q $path -i 99.5 -l 1000 > sample_nucmer/$file.filter
    show-coords -rclT $NUC_DIR/$file.filter > sample_nucmer/$file.coords
  done &
```

### Calculate the chain files:
```shell
  /scripts/crossmap_delta_to_chain.pl -fwd work/all_FWD.chain -rev work/all_REV.chain work/all_r.v1.q.v2.filter
```

### Use the chain files to map features in a supplied BED file. Do this for both the FWD and REV chain files.
```shell
  CrossMap.py bed work/all_FWD.chain sample_snps/SNP50K_Wm82.a1.bed | grep -v Fail > work/all_map_a1_a2_fwd.tsv
  CrossMap.py bed work/all_REV.chain sample_snps/SNP50K_Wm82.a1.bed | grep -v Fail > work/all_map_a1_a2_rev_TMP
```

### Make file of query chromosome & scaffold sizes (derived from filter files)
```shell
	grep -h '>' sample_nucmer/*filter | sed 's/>//' | awk -v OFS="\t" '{print $2, $4}' |
    sort -u > work/Q_chr_sizes.tsv
```

### Recover (translate) the reverse mappings
```shell
  scripts/crossmap_recover_rev_coords.pl -q work/Q_chr_sizes.tsv \
		-m work/all_map_a1_a2_rev_TMP -o work/all_map_a1_a2_rev.tsv
```

The final result (the feature mapping between genome assemblies) is the combination of the following two files: work/all_map_a1_a2_fwd.tsv work/all_map_a1_a2_rev.tsv

