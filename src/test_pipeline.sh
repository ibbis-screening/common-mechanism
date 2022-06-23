#! /usr/bin/env bash

# usage: ../../src/test_pipeline.sh # from inside path to folder of test sequences

# run the Common Mechanism on each test sequence

for i in $1; do ../../src/run_pipeline.sh $i; done

# collate test results

../../test_scripts/parse_test_res.py
