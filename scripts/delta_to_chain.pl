#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my $help;

GetOptions (
  "help:s" =>            \$help
);

my $usage = <<EOS;
  Usage:  delta_to_chain.pl *.filter

  Given one or more files in MUMMer delta file (.delta or .filter), generate
  a file in the UCSC "chainfile" format, for use by CrossMap.py

  -help:          display this help message
EOS

die "\n$usage\n" if ($help);
die "\n$usage\n!! Please provide one or more MUMMer delta or filter files via ARGV\n\n" if (@ARGV==0);


#chain printout components
my ($tName,$tSize,$tStrand,$tStart,$tEnd);
my ($qName,$qSize,$qStrand,$qStart,$qEnd);
my ($ca,$ct,$cq) = (0,0,0);
my $size = -1;
my ($dt,$dq) = (0,0);
my $score;
my $id=1;

foreach(@ARGV) {
  my $filename = $_;
  open (my $fh, '<', $filename) or die;

  while(<$fh>){
      my @line = split(/\s+/, $_);
      if($line[0] =~ /^>.*/){
        ## Delta header starting with > 
        ($tName, $qName, $tSize, $qSize) = @line;
        $tName =~s/^>//;
      } elsif (scalar @line > 2){
        ## set remainder of chain header using line following >
        my ($score2, $frame); # last elements; usless in this context
        ($tStart, $tEnd, $qStart, $qEnd, $score, $score2, $frame) = @line;
        $tStrand = $tStart < $tEnd ? '+' : '-';
        $qStrand = $qStart < $qEnd ? '+' : '-';
        if($qStrand eq '-'){ # if negative orientation, swap qStart and qEnd
          my $tmp = $qStart;
          $qStart = $qEnd;
          $qEnd = $tmp;
        }
        printf "%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\n",
          'chain',$score,$tName,$tSize,$tStrand,$tStart,$tEnd,$qName,$qSize,$qStrand,$qStart,$qEnd,$id;
        $id++;
      } else {
        if(!($line[0] =~ /^-?\d+\.?\d*$/)) {
          ## skip non-numeric lines (file headers, empty space);
          next;
        }
        my $dist = $line[0];
        if(abs($dist) > 1 || ($dist != 0 && $size == -1)){
          ## Set size
          if($size > -1){ printf("%d\t%d\t%d\n",$size,$dt,$dq);}
          $size = abs($dist)-1;
          $ca += $size;
          if($dist < 0){
            $dq = 1;
            $dt = 0;
            $cq++;
          } else {
            $dq = 0;
            $dt = 1;
            $ct++;
          }
        } elsif ( $dist == -1){
          $dq++;
          $cq++;
        } elsif ($dist == 1){
          $dt++;
          $ct++;
        } else {
          ## $line[0] = 0, end of current delta read print and reset counters;
          if($size > -1){printf("%d\t%d\t%d\n",$size,$dt,$dq);}
          printf("%d\n\n", abs($tEnd-($tStart+$ca+$ct)));
          $size = -1;
          ($ca,$cq,$ct) = (0,0,0);
          printf "\n";
      }
    }
  }
}

__END__

VERSIONS
AW = Andrew Wilkey; AB = Annie Brown; SC = Steven Cannon

v02 AW 2018-05-01 initial version 
v03 AB 2018-06-08 added a statement to insert a blank line at the end of each chain
v05 SC 2018-08-05 from deltaToChain-rev4-072518.pl to de-swap qEnd qStart in inversions
v06 SC 2018-08-06 refactor; add usage info

