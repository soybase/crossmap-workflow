# crossmap-workflow

## Basic workflow:
Run nucmer on each chromosome and scaffold against the whole genome assembly (to allow each region
to find its best home on the whole genome), then filter using the MUMmer utility delta-filter
(recommend parameters for a lenient mapping between genome assemblies: -i 99.5 -l 1000).

Then use the script delta_to_chain.pl in this repository to generate a "chain" file from the
delta-format .filter files.

Then use the CrossMap.py script (http://crossmap.sourceforge.net) to generate a mapping of a file of genomic coordinates. The instructions below assume local installation of the MUMmer package and CrossMap. Note that CrossMap can be easily installed using Conda, e.g.
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

### Calculate the chain file:
```shell
  ./scripts/delta_to_chain.pl sample_nucmer/*filter \
    > sample_out/nucmer_r.v1.q.v2.id995_len1k.chain
```

### Calculate the mapping with CrossMap
```shell
  CrossMap.py bed sample_out/nucmer_r.v1.q.v2.id995_len1k.chain sample_snps/soysnp50k.a1.pos_named \
    > sample_out/nucmer_r.v1.q.v2.id995_len1k.mapped_named
```


