#!/bin/bash

BAM=$1
SAMPLE=$2
BED=$3
N=$4
seq2cov.pl -b $BAM -N $SAMPLE $BED > cov.txt.${N}
touch seq2c.done.${N}
