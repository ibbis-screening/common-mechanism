from utils import *
import sys, os
import pandas as pd

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

file = sys.argv[1] + ".biorisk.hmmsearch"

# read in HMMER output and check for valid hits
if checkfile(file) == 1:
    hmmer = readhmmer(file)
    hmmer = hmmer[hmmer['E-value']<1e-10]
    if hmmer.shape[0] > 0:
        print("Biorisks: FLAG")
    else:
        print("Biorisks: PASS")
elif checkfile(file) == 2:
	print("Biorisks: PASS")
else:
	print("Biorisks: unexpected outcome")
