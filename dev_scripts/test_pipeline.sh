#! /usr/bin/env bash
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
# usage: ../src/test_pipeline.sh # from inside path to folder of test sequences

# run the Common Mechanism on each test sequence
for i in $1; do ../src/run_pipeline.sh $i > ${i//fasta/screen}; done

# collate test results
../flag.py .
