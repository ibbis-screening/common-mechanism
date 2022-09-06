#! /Users/wheelern/opt/miniconda3/bin/python

# parse all .screen.txt files in a directory and count stats on each flag and overall calls
# Usage: parse_test_res.py (in directory of choice)

import glob
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

# columns in the summary file, to be filled with the checks below
names = []
biorisk = []
reg_virus = []
reg_bact = []
benign = []



def check_flags(line, bin_list):
        if len(matching) > 0:
            if 'PASS' in matching[0]:
                bin_list.append("P")
            elif 'FLAG' in matching[0]:
                bin_list.append("F")
#                print(matching)
            else:
                bin_list.append(None)
        else:
            bin_list.append(None)
        
        return bin_list
    
for res in glob.glob('*.screen'):
#    print(res)
    names.append(res)
    with open(res) as f:
        lines = f.readlines()
        
        # biorisk screen
        matching = [s for s in lines if "Biorisk" in s]
#        print(matching)
        biorisk = check_flags(matching, biorisk)

        # reg_virus screen
        matching = [s for s in lines if "Regulated virus top hit: FLAG" in s]
#        print(matching)
        reg_virus = check_flags(matching, reg_virus)
        
        # reg_virus screen
        matching = [s for s in lines if "Regulated bacteria top hit: FLAG" in s]
#        print(matching)
        reg_bact = check_flags(matching, reg_bact)
        
        # benign screen - 1 means a regulated region failed to clear, 0 means benign coverage and clear
        nohits = [s for s in lines if "No housekeeping genes found" in s]
        fail = [s for s in lines if "Regulated region failed to clear" in s]
        toolow = [s for s in lines if "Housekeeping genes - <90% coverage achieved = FAIL" in s]
        passs = [s for s in lines if "Housekeeping genes - >90% coverage achieved" in s]
        clear = [s for s in lines if "No regulated regions to clear" in s]
        # if any region failed to clear, keep flag
        if len(nohits) > 0 or len(fail) > 0 or len(toolow) > 0:
            benign.append("F")
        # if none failed and passes are observed, drop flag
        elif len(passs) > 0:
            benign.append("P")
        elif len(clear) > 0:
            benign.append("-")
        else:
            benign.append(None)

#print(len(names), len(biorisk), len(reg_virus), len(reg_bact), len(benign))

breakdown = list(zip(names, biorisk, reg_virus, reg_bact, benign))
summary = []
for name, risk, reg_vir, reg_bac, ben in breakdown:
        # if a biorisk is flagged, flag the whole thing
    if risk == "F":
        summary.append((name, "F"))
#        print("Biorisk found")
    elif reg_vir == "F":
        summary.append((name, "F"))
#        print("Regulated virus found")
    # if it's a regulated bacterial pathogen but a known benign gene, clear it
    elif (reg_bac == "F" and ben == "P") == 1:
        summary.append((name, "P"))
#        print("Regulated bacterial housekeeping found")
    # if it's a regulated bacterial hit, flag it
    elif reg_bac == "F":
#        print("Nothing found")
        summary.append((name, "F"))
    else:
        summary.append((name, None))
#print(summary)
pd.DataFrame(summary).to_csv("test_summary.csv", index=False, header=None)

breakdown = pd.DataFrame(breakdown)
breakdown.columns = ("names", "biorisk", "regulated_virus", "regulated_bacteria", "benign")
breakdown.to_csv("test_itemized.csv", index=False)


#g = sb.FacetGrid(breakdown, col="virus")
#g.map_dataframe(sb.stripplot, x=breakdown["biorisk"], y=breakdown["reg_tax"], hue=breakdown["benign"])
#plt.savefig("Positive_set.png")

#sb.set_context("talk")
#plt.title("Predictions on test set")
##sb.stripplot(x=biorisk, y=reg_tax, hue=benign, data=breakdown, jitter=0.3, dodge=True, size=10)
#sb.swarmplot(x=biorisk, y=reg_bact, hue=benign, data=breakdown, dodge=True, size=10)
#sb.despine()
#plt.ylabel("Regulated pathogen hit")
#plt.xlabel("Biorisk hit")
#plt.xticks(rotation=30, ha='right')
#plt.savefig("Test_set.png",bbox_inches='tight')

print((pd.DataFrame(summary)[1]=="F").sum(), "/", len(summary))
