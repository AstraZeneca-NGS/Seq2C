#!/bin/bash
set -x

BAM=$1
SAMPLE=$2
BED=$3
N=$4
SEQ2COV=$5

SAMTOOLS='samtools'
if [ $6 ]
    then
    SAMTOOLS=$6
fi

OUTPUT=cov.txt.${N}
if [ $7 ]
    then
    OUTPUT=$7
fi

DONE_MARKER=seq2c.done.${N}
if [ $8 ]
    then
    DONE_MARKER=$8
fi

perl $SEQ2COV -m $SAMTOOLS -z -b $BAM -N $SAMPLE $BED > $OUTPUT
touch $DONE_MARKER

set +x