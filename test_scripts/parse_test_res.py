#! /Users/wheelern/opt/miniconda3/bin/python

# parse all .screen.txt files in a directory and count stats on each flag and overall calls

import glob
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

names = []
virus = []
biorisk = []
reg_tax = []
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
    
for res in glob.glob('*.screen.txt'):
#    print(res)
    names.append(res)
    with open(res) as f:
        lines = f.readlines()
        
        # biorisk screen
        matching = [s for s in lines if "Biorisk" in s]
#        print(matching)
        biorisk = check_flags(matching, biorisk)

        # reg_tax screen
        matching = [s for s in lines if "Regulated pathogens" in s]
#        print(matching)
        reg_tax = check_flags(matching, reg_tax)
        
        # is it a virus?
        matching = [s for s in lines if "virus" in s]
        if len(matching) > 0:
            virus.append("F")
        else:
            virus.append("P")
        
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

print(len(names), len(virus), len(biorisk), len(reg_tax), len(benign))

breakdown = list(zip(names, virus, biorisk, reg_tax, benign))
summary = []
for name, viral, risk, reg, ben in breakdown:
    if viral == "F" and reg == "F":
        summary.append((name, "F"))
    elif not None in (risk, reg, ben):
        if risk == "F":
            summary.append((name, "F"))
        elif (reg == "F" and ben == "P") == 1: # if it's a regulated bacterial pathogen but a known benign gene, clear it
            summary.append((name, "P"))
        elif reg == "F" and risk == "F":
            summary.append((name, "F"))
    else:
        summary.append((name, None))
pd.DataFrame(summary).to_csv("test_summary.csv", index=False, header=None)

breakdown = pd.DataFrame(breakdown)
breakdown.columns = ("names", "virus", "biorisk", "reg_tax", "benign")
breakdown.to_csv("test_itemized.csv", index=False)


#g = sb.FacetGrid(breakdown, col="virus")
#g.map_dataframe(sb.stripplot, x=breakdown["biorisk"], y=breakdown["reg_tax"], hue=breakdown["benign"])
#plt.savefig("Positive_set.png")

sb.set_context("talk")
plt.title("Predictions on test set")
#sb.stripplot(x=biorisk, y=reg_tax, hue=benign, data=breakdown, jitter=0.3, dodge=True, size=10)
sb.swarmplot(x=biorisk, y=reg_tax, hue=benign, data=breakdown, dodge=True, size=10)
sb.despine()
plt.ylabel("Regulated pathogen hit")
plt.xlabel("Biorisk hit")
plt.xticks(rotation=30, ha='right')
plt.savefig("Test_set.png",bbox_inches='tight')

print((pd.DataFrame(summary)[1]=="F").sum(), "/", len(summary))
