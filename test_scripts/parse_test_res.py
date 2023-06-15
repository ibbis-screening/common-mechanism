#! /usr/bin/env python

# parse all .screen.txt files in a directory and count stats on each flag and overall calls
# Usage: parse_test_res.py (in directory of choice)

import glob
import pandas as pd

# columns in the summary file, to be filled with the checks below
names = []
biorisk = []
vfs = []
reg_virus = []
reg_bact = []
reg_fungi = []
reg_nonreg = []
benign = []

def check_flags(matching, bin_list):
    if len(matching) > 0:
        hit = 0
        for match in matching:
            if hit == 0:
                if ': FLAG' in match:
                    hit = 1
                    bin_list.append("F")
        if hit == 0:
            bin_list.append("P")
            # print(matching)
    else:
        bin_list.append("Err")
    
    return bin_list
    
for res in glob.glob('*.screen'):
    names.append(res)
    # print(res)
    with open(res) as f:
        lines = f.readlines()
        
        # biorisk screen
        matching = [s for s in lines if "Regulated gene" in s]
        biorisk = check_flags(matching, biorisk)

        # VF screen
        vf = [s for s in lines if "Virulence factor found" in s]
        if len(vf) > 0:
            vfs.append("F")
        else:
            vfs.append("P")
        
        # homol_fail = [s for s in lines if "Homology search has failed" in s]
        # if len(homol_fail) > 0:
        #     reg_virus.append("Err")
        #     reg_bact.append("Err")
        #     reg_fungi.append("Err")
        # else:

        # reg_virus screen - fetch all coding and noncoding reports
        matching_virus = [s for s in lines if "found in only regulated organisms: FLAG (virus)" in s]
        # print(matching_virus)
        if len(matching_virus) > 0:
            reg_virus = check_flags(matching_virus, reg_virus)
        else:
            reg_virus.append("P")
        
        # reg_bact screen - fetch all coding and noncoding reports
        matching_bact = [s for s in lines if "found in only regulated organisms: FLAG (bacteria)" in s]
        # print(matching_bact)
        if len(matching_bact) > 0:
            reg_bact.append("F")
            # reg_bact = check_flags(matching_bact, reg_bact)
        else:
            reg_bact.append("P")

        # reg_fungi screen - fetch all coding and noncoding reports
        matching_fungi = [s for s in lines if "found in only regulated organisms: FLAG (fungi)" in s]
        # print(matching_fungi)
        if len(matching_fungi) > 0:
            reg_fungi = check_flags(matching_fungi, reg_fungi)
        else:
            reg_fungi.append("P")

        # reg_nonreg screen - ID cases where the same sequence is found in regulated and non-regulated organisms
        matching_reg_nonreg = [s for s in lines if "found in both regulated and non-regulated organisms" in s]
        if len(matching_reg_nonreg) > 0:
            reg_nonreg.append("F")
        else:
            reg_nonreg.append("P")

        # benign screen - 1 means a regulated region failed to clear, 0 means benign coverage and clear
        allpass = [s for s in lines if "all regulated regions cleared: PASS" in s]
        anyfail = [s for s in lines if "failed to clear: FLAG" in s]
        clear = [s for s in lines if "no regulated regions to clear" in s]
        # if any region failed to clear, keep flag
        if len(allpass) > 0:
            benign.append("P")
        # if none failed and passes are observed, drop flag
        elif len(anyfail) > 0:
            benign.append("F")
        elif len(clear) > 0:
            benign.append("-")
        else:
            benign.append("Err")
    # print(len(biorisk), len(reg_virus), len(reg_bact), len(benign))

#print(len(names), len(biorisk), len(reg_virus), len(reg_bact), len(benign))

breakdown = list(zip(names, biorisk, vfs, reg_virus, reg_bact, reg_fungi, reg_nonreg, benign))
summary = []
for name, risk, vf, reg_vir, reg_bac, reg_fungi, reg_nonreg, ben in breakdown:
        # if a biorisk is flagged, flag the whole thing
    if risk == "Err" or reg_bac == "Err" or ben == "Err":
        summary.append((name, "Err"))
    else:
        if risk == "F":
            summary.append((name, "F"))
    #        print("Biorisk found")
        elif reg_vir == "F":
            summary.append((name, "F"))
    #        print("Regulated virus found")
        # if it's a regulated bacterial pathogen but a known benign gene, clear it
        elif (reg_bac == "F" and ben == "P" and vf == "P") == 1:
            summary.append((name, "P"))
        elif (reg_bac == "F" and ben == "P" and vf == "F") == 1: # some VFs are also housekeeping genes, but these seem to be poorly supported Victors genes
            summary.append((name, "P"))
    #        print("Regulated bacterial housekeeping found")
        elif (reg_fungi == "F" and ben == "P") == 1:
            summary.append((name, "P"))
        # if it's a regulated bacterial hit, flag it
        elif reg_bac == "F":
            summary.append((name, "F"))
        elif reg_fungi == "F":
            summary.append((name, "F"))
        else:
            summary.append((name, "P"))
pd.DataFrame(summary).to_csv("test_summary.csv", index=False, header=None)

breakdown = pd.DataFrame(breakdown)
breakdown.columns = ("names", "biorisk", "virulence factor", "regulated_virus", "regulated_bacteria", "regulated_fungi/oomycetes", "mix of reg and non-reg", "benign")
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

print("Flags: ", (pd.DataFrame(summary)[1]=="F").sum(), "/", len(summary))
print("Errors: ", (pd.DataFrame(summary)[1]=="Err").sum())
