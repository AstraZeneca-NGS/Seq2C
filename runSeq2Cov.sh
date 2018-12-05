#!/bin/bash

# Given a bam and bed file, calculate the coverage
# The script is used to submit jobs to grid

BAM=$1
SAMPLE=$2
BED=$3

~/work/NGS/NGS/seq2c/seq2cov.pl -N $SAMPLE -b $BAM $BED > $SAMPLE.cov
