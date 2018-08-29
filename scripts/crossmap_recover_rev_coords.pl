#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my ($q_sizes, $map, $out, $help);

GetOptions ( 
  "q_sizes=s" => \$q_sizes,
  "map=s"     => \$map,
  "out=s"     => \$out,
  "help"      => \$help 
);

my $usage = <<EOS;
  Usage: crossmap_recover_rev_coords.pl -q_sizes FILE -map FILE

  Given a file of query chromosome (and scaffold) sizes and a file of position mappings
  from CrossMap, for inverted regions (output from crossmap_delta_to_chain.pl), 
  generate a mapping with coordinates transformed as: (q_chromosome minus new_position)

  The reason for this step (and script): as of late 2018, it seems that liftOver and CrossMap 
  are doing the wrong thing with inversions. The script crossmap_delta_to_chain.pl generates 
  chain-format output from MUMmer input, separated into forward and reverse chain files
  (one for forward alignments and one for inversions). Then CrossMap is run separately on both,
  and then this script recovers the correct forward mapping for the features that fall in
  inversions (via the reverse mappings).

  Required:
  -q_sizes  (required) filename for chromosome/scaffold and sizes (two columns)
  -map  (required) filename for CrossMap output, for inverted regions, after
              running crossmap_delta_to_chain.pl
  -out      output of the mapping translation, in CrossMap output format

  Options:
  -help:    display this help message
EOS

die "\n$usage\n" if ($help);
die "\n$usage\n" unless ($q_sizes && $map && $out);

open my $QSZ, "<", $q_sizes or die "can't open in $q_sizes: $!\n";
open my $MAP, "<", $map or die "can't open in $map: $!\n";
open my $OUT, ">", $out or die "can't open out $out: $!\n";

my %Q_sizes;
while (<$QSZ>) {
  chomp;
  next if $_ =~ /^$/;
  my ($chr, $size) = split(/\s/, $_);
  #print "$chr\t$size\n";
  $Q_sizes{$chr} = $size;
}

while (<$MAP>) {
  chomp;
  next if $_ =~ /^$/;
  my ($chrR, $startR, $endR, $idR, $arrow, $chrQ, $startQ, $endQ, $idQ) = split(/\s/, $_);

#print "==", $Q_sizes{$chrQ}, "\t$chrQ\t$startQ\t$endQ\n";

  my $new_startQ = $Q_sizes{$chrQ} - $startQ;
  my $new_endQ   = $Q_sizes{$chrQ} - $endQ;
  print $OUT "$chrR\t$startR\t$endR\t$idR\t$arrow\t$chrQ\t$new_startQ\t$new_endQ\t$idQ\n";
}

__END__

VERSIONS
SC=Steven Cannon

2018:
v01 08-22 SC initial version
v02 08-29 SC fix hashing - switching from $chrR to $chrQ


