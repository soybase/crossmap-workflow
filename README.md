# crossmap-workflow

## Basic workflow:
Run nucmer on each chromosome and scaffold against the whole genome assembly (to allow each region
to find its best home on the whole genome), then filter using the MUMmer utility delta-filter
(recommend parameters for a lenient mapping between genome assemblies: -i 99.5 -l 1000).

Then use the script delta_to_chain.pl in this repository to generate a "chain" file from the
delta-format .filter files.

Then use the CrossMap.py script to generate a mapping of a file of genomic coordinates.

## Example:

### Filter the raw nucmer output:
  NUC_DIR="nucmer_output"
  for path in $NUC_DIR/*.delta; do
    file=`basename $path .delta`
    delta-filter -q $path -i 99.5 -l 1000 > $NUC_DIR/$file.filter
    show-coords -rclT $NUC_DIR/$file.filter > $NUC_DIR/$file.coords
  done &

### Calculate the chain file:
```shell
  ./scripts/delta_to_chain.pl nucmer_r.v1.q.v2.id995_len1k/*filter > nucmer_r.v1.q.v2.id995_len1k.chain
  wc nucmer_r.v1.q.v2.id995_len1k.chain
```

### Calculate the mapping with CrossMap
```shell
  CrossMap.py bed nucmer_r.v1.q.v2.id995_len1k.chain soysnp50k.a1.pos_named \
    > nucmer_r.v1.q.v2.id995_len1k.mapped_named
```

