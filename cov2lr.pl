#!/usr/bin/env perl

# Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium 
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).

use warnings;
use FindBin;
use lib "$FindBin::Bin/libraries/";
use Stat::Basic;
use Statistics::TTest;
use Getopt::Std;
use strict;

our ($opt_c, $opt_a, $opt_H, $opt_F);

getopts( 'Hac:F:' ) || USAGE();
$opt_H && USAGE();
my $FAILEDFACTOR = $opt_F ? $opt_F : 0.2;

my $CNT = shift;
my %cnt;
open( CNT, $CNT );
my @cnt;
while(<CNT>) {
    chomp;
    next if ( /Undetermined/ );
    next if ( /Sample\t/ || /Total/ );
    my @a = split(/\t/);
    $cnt{$a[0]} = $a[1];
    push(@cnt, $a[1]);
}
close(CNT);
my $stat = new Stat::Basic;
my $meanreads = $stat->prctile(\@cnt, 50);
my %factor; # to adjust sequencing coverage
while(my ($k, $v) = each %cnt) {
    $factor{ $k } = $meanreads/$v;
}
my %raw; # key: gene value: hash of (key: sample; value: raw depth)
my %norm1; # key: gene value: hash of (key: sample; value: normalized depth by sequencing distribution)
my %norm1b; # key: gene value: hash of (key: sample; value: normalized depth by gene)
my %norm2; # key: gene value: hash of (key: sample; value: normalized by gene median scaling log2)
my %norm3; # key: gene value: hash of (key: sample; value: normalized by sample median scaling log2)
my %norm_c;
my %samples;
my @depth;
my %data;
my %loc;
while( <> ) {
    next if ( /Whole-/ );
    next if ( /_CONTROL_/ );
    next if ( /^Sample\t/ );
    my ($sample, $gene, $chr, $s, $e, $desc, $len, $depth) = split(/\t/);
    next if ( $sample eq "Undetermined" );
    next unless( $cnt{ $sample } );
    my $k = $opt_a ? join("\t", $gene, $chr, $s, $e, $len) : $gene;
    $loc{ $gene }->{ chr } = $chr;
    $loc{ $gene }->{ start } = $s unless( $loc{ $gene }->{ start } && $loc{ $gene }->{ start } < $s );
    $loc{ $gene }->{ end } = $e unless( $loc{ $gene }->{ end } && $loc{ $gene }->{ end } > $e );
    $loc{ $gene }->{ len } += $len;
    $data{ $k }->{ $sample }->{ len } += $len;
    $data{ $k }->{ $sample }->{ cov } += $depth * $len;
}

while(my($k, $v) = each %data) {
    while( my($sample, $dv) = each %$v ) {
	$raw{ $k }->{ $sample } = $dv->{cov}/$dv->{len};
	$norm1{ $k }->{ $sample } = sprintf("%.2f", $raw{ $k }->{ $sample }*$factor{ $sample });
	$samples{ $sample } = 1;
	push(@depth, $norm1{ $k }->{ $sample } );
    }
}

my @samples = keys %samples;
my $meddepth = $stat->median(\@depth);

# remove genes/amplicons that failed
my %bad;
my @gooddepth = ();
while(my($k, $v) = each %data) {
    my @tmp = map { $norm1{ $k }->{ $_} } @samples;
    my $kp80 = $stat->prctile(\@tmp, 80);
    if ( $kp80 < $meddepth*$FAILEDFACTOR ) {
        $bad{ $k } = 1;
    } else {
        push(@gooddepth, @tmp);
    }
}
$meddepth = $stat->median(\@gooddepth); # re-adjust median depth using only those from good amplicons/genes
#$meddepth = $stat->prctile(\@gooddepth, 25); # re-adjust median depth using only those from good amplicons/genes

my %factor2;  # Gene/amplicon factor
while( my ($k, $v) = each %norm1) {
    next if ( $bad{ $k } );
    my @t = values %$v;
    my $smt = $stat->median(\@t);
    $factor2{ $k } = $smt != 0 ? $meddepth/$smt : 0;
}

while( my ($k, $v) = each %norm1) {
    next if ( $bad{ $k } );
    foreach my $s (@samples) {
	$norm1b{ $k }->{ $s } = sprintf("%.2f", $v->{$s} * $factor2{ $k }+0.1);
        $norm2{ $k }->{ $s } = sprintf("%.2f", log(($v->{$s} * $factor2{ $k }+0.1)/$meddepth)/log(2));
        #$norm3{ $k }->{ $s } = sprintf("%.2f", log(($v->{$s} * $factor2{ $k }+0.1)/$samplemedian{ $s })/log(2));
    }
}

my %samplemedian;
foreach my $s (@samples) {
    my @tmp = ();
    while( my ($k, $v) = each %norm1b ) {
	next if ( $bad{ $k } ); # ignore failed genes/amplicons
        push( @tmp, $v->{ $s } );
    }
    $samplemedian{ $s } = $stat->median( \@tmp );
}

while( my ($k, $v) = each %norm1b) {
    next if ( $bad{ $k } );
    foreach my $s (@samples) {
        $norm3{ $k }->{ $s } = sprintf("%.2f", log($v->{$s}/$samplemedian{ $s })/log(2));
    }
}

my %factor_c;  # Gene/amplicon factor
if ( $opt_c ) {
    while( my ($k, $v) = each %norm1b) {
        next if ($bad{ $k });
        while( my ($s, $d) = each %$v ) {
            my @controls = split(/:/, $opt_c);
            my @tcntl = map { $norm1b{ $k }->{ $_ } } @controls;
            my $cntl_ave = $stat->mean( \@tcntl );
            $norm_c{ $k }->{ $s } = sprintf("%.3f", log($norm1b{ $k }->{ $s }/$cntl_ave)/log(2));
        }
    }
    foreach my $s (@samples) {
        my @tmp = ();
        while( my ($k, $v) = each %norm_c ) {
            next if ( $bad{ $k } ); # ignore failed genes/amplicons
            push( @tmp, $v->{ $s } );
        }
        $factor_c{ $s } = $stat->median( \@tmp );
    }
}

while( my ($k, $v) = each %norm2) {
    next if ($bad{ $k });
    while( my ($s, $d) = each %$v ) {
	my $t = $opt_a ? $k : "$k\t$loc{$k}->{chr}\t$loc{$k}->{start}\t$loc{$k}->{end}\t$data{$k}->{$s}->{len}";
        my $str = join("\t", $s, $t, $raw{ $k }->{ $s }, $norm1{ $k }->{ $s }, $norm1b{ $k }->{ $s }, $d, $norm3{ $k }->{ $s });
	if ( $opt_c ) {
	    $str .= "\t" .  sprintf("%.3f", $norm_c{ $k }->{ $s } - $factor_c{ $s });
	}
	print "$str\n";
    }
}

exit(0);

sub USAGE {
getopts( 'aPc:F:s:A:D:' );
print <<USAGE;
    Usage: $0 [-aH] [-c control] mapping_reads coverage.txt

    The $0 program will convert a coverage file to copy number profile.

    Arguments are:
    mapping_reads: Required.  A file containing # of mapped or sequenced reads for samples.  At least two columns.
                   First is the sample name, 2nd is the number of mapped or sequenced reads.
    coverage.txt:  The coverage output file from checkCov.pl script.  Can also take from standard in or more than
                   one file.

    Options are:

    -a Indicate this is amplicon or exon based calling.  By default, it'll aggregate at gene level.

    -c sample(s)
       Specify the control sample(s), if aplicable.  Multiple controls are allowed, which are separated by ":"

    -F double
       The failed factor for individual amplicons.  If (the 80th percentile of an amplicon depth)/(the global median depth)
       is less than the argument, the amplicon is considered failed and won't be used in calculation.  Default: 0.2.

AUTHOR
       Written by Zhongwu Lai, AstraZeneca, Boston, USA

REPORTING BUGS
       Report bugs to zhongwu\@yahoo.com

COPYRIGHT
       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.

USAGE
exit(0);
}
