#!/bin/bash

BAM=$1
SAMPLE=$2
BED=$3
N=$4
SEQ2COV=$5
perl $SEQ2COV -z -b $BAM -N $SAMPLE $BED > cov.txt.${N}
touch seq2c.done.${N}
