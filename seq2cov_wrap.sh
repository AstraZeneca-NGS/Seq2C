#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
BAM=$1
SAMPLE=$2
BED=$3
N=$4
perl ${DIR}/seq2cov.pl -z -b $BAM -N $SAMPLE $BED > cov.txt.${N}
touch seq2c.done.${N}
