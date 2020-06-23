#!/usr/bin/env perl

use Getopt::Std;
use strict;

our ($opt_n, $opt_d, $opt_a, $opt_e, $opt_A, $opt_D, $opt_N, $opt_p, $opt_g, $opt_P, $opt_H, $opt_G, $opt_k);

getopts('Hkgn:d:a:e:D:A:N:p:G:P:') || USAGE();
$opt_H && USAGE();

my %purity;
if ( $opt_p ) {
    setPurity($opt_p);
}

my $MINDEL = $opt_D ? $opt_D : -2.0;
my $MINAMP = $opt_A ? $opt_A : 1.45;  # 6 copies for pure sample

my $MINEXONDEL = $opt_d ? $opt_d : -2.0;
my $MINEXONAMP = $opt_a ? $opt_a : 1.75;
my $MINEXON = $opt_e ? $opt_e : 1;
my $N = $opt_N ? $opt_N : 5; # If a breakpoint is called more than $N samples, then it's deemed a false positive and filtered
my @GAIN = $opt_G ? split(/:/, $opt_G) : ("MYC");
my %GAIN = map { ($_, 1); } @GAIN;

my %count;
my %samples;
if ( $opt_k ) {
    print join("\t", qw(Sample Empty Variant_Type Gene NA - - Segment - - Transf_LogRatio Segments LogRatio alteration - - - - - - - Alteration));
    print "\n";
}
while( <> ) {
    my @a = split(/\t/);
    my $gene = $a[1];
    my $sample = $a[0];
    if ( $opt_n ) {
	if ( $a[0] =~ /$opt_n/ ) {
	    $sample = $1; 
	} else {
	    next;
	}
    }
    $samples{ $sample } = 1;
    my ($SMINDEL, $SMINAMP, $SMINEXONDEL, $SMINEXONAMP) = ($MINDEL, $MINAMP, $MINEXONDEL, $MINEXONAMP);
    if ( $purity{ $sample } ) {
    } elsif ( $opt_P ) {
    }
    my $lr = $a[6];
    #if ( $a[10] && $a[10] < $MINEXON ) {
	#$lr = $a[12];
        #next unless( $lr > $MINAMP || $lr < $MINDEL );
    #}
    my $desc = ($a[8] && $a[8] eq "BP") ? "$a[10] of $a[11]" : "$a[11] of $a[11]";
    if ( $opt_g || $GAIN{ $gene } ) { # Only do whole gene for copy gains
        if ( $lr >= 0.75 && $lr < $SMINAMP ) {
	    print join("\t", $sample, "", "copy-number-alteration", $a[1], "NA", "-", "-", "$a[2]:$a[3]", "-", "-", sprintf("%.1f", 2**$lr*2), $desc, $lr, "gain", "-", "-", "-", "-", "-", "-", "-", "Gain"), "\n";
	    next;
	}
    }
    if ( $a[10] && $a[8] eq "BP" && ($a[10] >= $MINEXON || $a[11] - $a[10] >= $MINEXON) ) {
	    $lr = $a[12];
        next unless( $lr >= $SMINEXONAMP || $lr <= $SMINEXONDEL );
    }
    next unless( $lr >= $SMINAMP || $lr <= $SMINDEL );
    my $type = $lr <= $SMINDEL ? "Deletion" : "Amplification";  # Fixed bug when logR=-2.00
    next if ( $a[15] && $a[15] >= $N );
    print join("\t", $sample, "", "copy-number-alteration", $a[1], "NA", "-", "-", "$a[2]:$a[3]", "-", "-", sprintf("%.1f", 2**$lr*2), $desc, $lr, $type eq "Deletion" ? "loss" : "amplification", "-", "-", "-", "-", "-", "-", "-", $type), "\n";
}

# Set the tumor purity for each sample
sub setPurity {
    my $file = shift;
    open(PUR, $file);
    while( <PUR> ) {
        chomp;
	my ($sample, $purity) = split(/\t/);
	$purity /= 100 if ( $purity >= 5 );
	$purity{ $sample } = $purity;
    }
    close(PUR);
}

sub USAGE {
    print STDERR<<USAGE;
    $0 [-k] [-g] [-e exons] [-n reg] [-N num] [-A amp] [-a amp] [-D del] [-d del] [-p purity_file] [-P purity] lr2gene_output

    The program will parse seq2c output and make calls for each gene and output in the format compatible with OncoPrint.
    It has the option to provide the purity so that log2ratio thresholds will be adjusted accordingly.  By default, it calls
    genes that are homozygously deleted or amplified with >= 6 copies.

    Options:
    -k Print header

    -g Whether to output copy gains [4-5] copies.  Default: no.

    -p file
       A file contains the tumor purity for all samples.  Two columns, first is sample name, second is the purity in % or fraction [0-1].

    -P double
       The purity.  Default: 1 or 100%, as is for cell lines.  If set, all samples will assume to have the same purity.

    -n regex
       The regular expression to extract sample names.  Default: none.

    -N num
       If an breakpoint has been called in >= num of samples, it's deemed false positive.  Default: 5

    -e exons
       Minimum number of exons/amplicon.  Default: 1

    For whole gene:
    -D log2ratio
       The log2ratio to determine that a gene is homozygously deleted.  Default: -2.0

    -A log2ratio
       The log2ratio to determine that a gene is amplified.  Default: 1.45.

    For < 3 exons:
    -d log2ratio
       The minimum log2ratio to determine that 1-2 exons are deleted. Should be lower than [-d] to reduce false positives. Default: -2.5

    -a log2ratio
       The minimum log2ratio to determine that 1-2 exons are amplified. Should be larger than [-a] to reduce false positives. Default: 1.75
       
    For gains:
    -G Genes
       List of genes, seperated by ":", for which gain will be captured.  Default: MYC

AUTHOR
    Written by Zhongwu Lai

USAGE
    exit(0);
}
