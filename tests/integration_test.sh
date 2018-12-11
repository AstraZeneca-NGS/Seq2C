#!/bin/bash
OUTPUT="seq2c_results.txt"

#Run samples
echo Run sample1 test without control sample
SAMPLE2BAM="samples/Colo_sample1_hg19_chr20_sample2bam.txt"
BED="bed/hg19_chr20_Colo_sample1.bed"

seq2c.sh $SAMPLE2BAM $BED
ACTUAL_SAMPLE1="seq2c_results_sample1.txt"
mv $OUTPUT $ACTUAL_SAMPLE1


echo Run sample1 test with control sample
CONTROLSAMPLE="Colo829.chr20.18.Sample1"
SAMPLE2BAM="samples/Colo_sample1_hg19_chr20_sample2bam.txt"
BED="bed/hg19_chr20_Colo_sample1.bed"

seq2c.sh $SAMPLE2BAM $BED $CONTROLSAMPLE
ACTUAL_SAMPLE1_CONTROL="seq2c_results_sample1_control.txt"
mv $OUTPUT $ACTUAL_SAMPLE1_CONTROL

#Compare results
echo Compare results of sample 1
DIFF_FILE_SAMPLE1="diff_sample1.txt"
EXPECTED_SAMPLE1="expected/expected_Colo_sample1.txt"

if diff <(sort $EXPECTED_SAMPLE1) <(sort $ACTUAL_SAMPLE1) > $DIFF_FILE_SAMPLE1
then
	echo "OK: no differences";
else
	echo "ERROR: files have differences!"
	exit 1;
fi

echo Compare results of sample 1 with control sample
DIFF_FILE_SAMPLE1_CONTROL="diff_sample1_control.txt"
EXPECTED_SAMPLE1_CONTROL="expected/expected_Colo_sample1_control.txt"

if diff <(sort $EXPECTED_SAMPLE1_CONTROL) <(sort $ACTUAL_SAMPLE1_CONTROL) > $DIFF_FILE_SAMPLE1_CONTROL
then
	echo "OK: no differences";
else
	echo "ERROR: files have differences!"
	exit 1;
fi