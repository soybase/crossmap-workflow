#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my ($fwd_out, $rev_out, $help);

GetOptions ( 
  "fwd_out=s" => \$fwd_out,
  "rev_out=s" => \$rev_out,
  "help"      => \$help 
);

my $usage = <<EOS;
  Usage: crossmap_delta_to_chain.pl -fwd FILE -rev FILE *.filter

  Given one or more files in MUMMer delta file (.delta or .filter), generate
  a mapping in the UCSC "chain" format, for use by CrossMap.py

  As of late 2018, it seems that liftOver and CrossMap are doing the wrong thing
  with inversions, so this script produces two chain files: one for forward
  alignments and one for inversions. The one for inversions transforms the 
  coords as the Q chromosome length minus the original coordinates, in + orientation.
  The results need to be recovered afterwards as Qchr size - transated coordinates.

  Required:
  -fwd_out    (required) filename for forward chain components
  -rev_out    (required) filename for reverse chain components

  Options:
  -help:      display this help message
EOS

die "\n$usage\n" if ($help);
die "\n$usage\n!! Please provide one or more MUMMer delta or filter files via ARGV\n\n" if (@ARGV==0);

open my $FWD_OUT, ">", $fwd_out or die "can't open out $fwd_out: $!\n";
open my $REV_OUT, ">", $rev_out or die "can't open out $rev_out: $!\n";

#chain output components
my ($tName,$tSize,$tStrand,$tStart,$tEnd);
my ($qName,$qSize,$qStrand,$qStart,$qEnd);
my ($sum_aligns,$sumT,$sumQ) = (0,0,0);
my $size = -1;
my ($dt,$dq) = (0,0);
my $score;
my $id=1;

foreach(@ARGV) {
  my $filename = $_;
  open (my $fh, '<', $filename) or die "can't open in $filename: $!\n";

  while(<$fh>) {
    next if ( $_ =~ /^$/ );
    my @line = split(/\s+/, $_);
    if($line[0] =~ /^>.*/){ # Delta header starting with > 
      ($tName, $qName, $tSize, $qSize) = @line;
      $tName =~s/^>//;
    } 
    elsif (scalar @line == 7){ # set remainder of chain header using line following >
      my ($score2, $frame); # last elements; usless in this context
      ($tStart, $tEnd, $qStart, $qEnd, $score, $score2, $frame) = @line;
      $tStrand = $tStart < $tEnd ? '+' : '-';
      $qStrand = $qStart < $qEnd ? '+' : '-';
      if ($qStrand eq '-') { # in inversion; Swap $qStart, $qEnd
        # NOTE: translate query coords: qSize-coord, and change orientation to +
        printf $REV_OUT "%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n",
          'chain',$score,$tName,$tSize,$tStrand,$tStart,$tEnd,
                         $qName,$qSize,"+",$qSize-$qStart,$qSize-$qEnd,$id;
          #'chain',$score,$tName,$tSize,$tStrand,$tStart,$tEnd,
          #               $qName,$qSize,$qStrand,$qEnd,$qStart,$id;
      }
      else { # not in inversion; print $qStart, $qEnd
        printf $FWD_OUT "%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n",
          'chain',$score,$tName,$tSize,$tStrand,$tStart,$tEnd,
                         $qName,$qSize,$qStrand,$qStart,$qEnd,$id;
      }
      $id++;
    } 
    else { # delta section; two or one fields, or blank
      if(!($line[0] =~ /^-?\d+\.?\d*$/)) { # skip non-numeric lines (file headers, empty space);
        next;
      }
      my $dist = $line[0];
      if(abs($dist) > 1 || ($dist != 0 && $size == -1)) { 
        if ($qStrand eq '-') { # in inversion
          if($size > -1) { printf $REV_OUT "%d\t%d\t%d\n", $size,$dt,$dq }
        }
        else { # not in inversion
          if($size > -1) { printf $FWD_OUT "%d\t%d\t%d\n", $size,$dt,$dq }
        }

        $size = abs($dist)-1;
        $sum_aligns += $size;
        if($dist < 0) {
          $dq = 1;
          $dt = 0;
          $sumQ++;
        } 
        else {
          $dq = 0;
          $dt = 1;
          $sumT++;
        }
      } 
      elsif ( $dist == -1) {
        $dq++;
        $sumQ++;
      } 
      elsif ($dist == 1) {
        $dt++;
        $sumT++;
      } 
      else { # $line[0] = 0, end of current delta read print and reset counters;
        if ($qStrand eq '-') { # in inversion
          if ($size > -1) {
            printf $REV_OUT "%d\t%d\t%d\n", $size,$dt,$dq;
          }
        }
        else { # not in inversion
          if ($size > -1) {
            printf $FWD_OUT "%d\t%d\t%d\n", $size,$dt,$dq;
          }
        }

        if ($qStrand eq '-') { # in inversion
          printf $REV_OUT "%d\n\n", $qStart-($qEnd+$sum_aligns+$sumQ);
        }
        else { # not in inversion
          printf $FWD_OUT "%d\n\n", $qEnd-($qStart+$sum_aligns+$sumQ);
        }

        $size = -1;
        ($sum_aligns,$sumQ,$sumT) = (0,0,0);
      }
    }
  }
}

__END__

VERSIONS
AW = Andrew Wilkey; AB = Annie Brown; SC = Steven Cannon

2018:
v02 05-01 AW initial version 
v03 06-08 AB added a statement to insert a blank line at the end of each chain
v05 08-05 SC from deltaToChain-rev4-072518.pl to de-swap qEnd qStart in inversions
v06 08-15 SC Add usage info, and refactor and do significant rewrite. 
v07 08-15 SC Print forward and reverse components to separate chain files

