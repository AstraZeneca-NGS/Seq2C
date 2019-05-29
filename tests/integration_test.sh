#!/bin/bash -eu
set -o pipefail

#-------------------------------------------------------------------------------
# IMPORTANT: Before run you have to set paths to the BAM files of 
#            GTL_16_5 and controls in samples/GTL_16_5_sample2bam.txt
#-------------------------------------------------------------------------------

OUTPUT="seq2c_results.txt"
OUTPUT_DIR="/tmp/seq2c_test"

SCRIPT_DIR=`cd "$(dirname $0)"; pwd -P`
SEQ2C_PATH="$SCRIPT_DIR/.."
TEST_DIR="$SCRIPT_DIR"
echo "SCRIPT_DIR=$SCRIPT_DIR"
echo "SEQ2C_PATH=$SEQ2C_PATH"

# Add seq2c to path
PATH="$PATH:$SEQ2C_PATH"

#---
# Compare results
# Arguments: test_name expected_resutls actual_results
#---
function diff_results {
	TEST_NAME="$1"
	EXPECTED="$2"
	RESULTS="$3"

	DIFF_FILE="$OUTPUT_DIR/test_$TEST_NAME.diff.txt"
	echo "Test $TEST_NAME: Compare results, creating diff file '$DIFF_FILE'"

	if diff <(sort $EXPECTED) <(sort $RESULTS) > $DIFF_FILE
	then
		echo "Test $TEST_NAME: OK"
	else
		echo "Test $TEST_NAME: ERROR"
		exit 1;
	fi
}

#---
# Test 1
#---
function test_1 {
	echo "Test 1: Run sample1 test without control sample"
	SAMPLE2BAM="$TEST_DIR/samples/Colo_sample1_hg19_chr20_sample2bam.txt"
	BED="$TEST_DIR/bed/hg19_chr20_Colo_sample1.bed"
	ACTUAL="$OUTPUT_DIR/test_1.seq2c_results_sample1.txt"
	EXPECTED="$TEST_DIR/expected/expected_Colo_sample1.txt"

	seq2c.sh $SAMPLE2BAM $BED
	mv $OUTPUT $ACTUAL
	diff_results 1 "$EXPECTED" "$ACTUAL"
}

#---
# Test 2
#---
function test_2 {
	echo "Test 2: Run sample1 test with control sample"
	CONTROL_SAMPLE="Colo829.chr20.18.Sample1"
	SAMPLE2BAM="$TEST_DIR/samples/Colo_sample1_hg19_chr20_sample2bam.txt"
	BED="$TEST_DIR/bed/hg19_chr20_Colo_sample1.bed"
	ACTUAL="$OUTPUT_DIR/test_2.seq2c_results_sample1_control.txt"
	EXPECTED="$TEST_DIR/expected/expected_Colo_sample1_control.txt"

	seq2c.sh $SAMPLE2BAM $BED $CONTROL_SAMPLE
	mv $OUTPUT $ACTUAL
	diff_results 2 "$EXPECTED" "$ACTUAL"
}

#---
# Test 3
#---
function test_3 {
	echo "Test 3: Run GTL_16_5 test with control sample"
	CONTROL_SAMPLE="control1:control2:control3:control4:control5:control6:"
	SAMPLE2BAM="$TEST_DIR/samples/GTL_16_5_sample2bam.txt"
	BED="$TEST_DIR/bed/panel_az600_chr7.bed"
	ACTUAL="$OUTPUT_DIR/test_3.seq2c_results_GTL_16_5.txt"
	EXPECTED="$TEST_DIR/expected/expected_GTL_16_5.txt"

	seq2c.sh $SAMPLE2BAM $BED $CONTROL_SAMPLE
	mv $OUTPUT $ACTUAL
	diff_results 3 "$EXPECTED" "$ACTUAL"
}

#---
# Main: Run all tests
#---

mkdir "$OUTPUT_DIR" || true

test_1
test_2
test_3

