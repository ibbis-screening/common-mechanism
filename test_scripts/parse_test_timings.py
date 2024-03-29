#! /usr/bin/env python3

# parse all .stats files in a folder to plot the distribution of time taken to screen
# Usage: parse_test_timings.py (in directory of choice)

# make it compatible with both individual files and a pooled one

import glob
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

from datetime import datetime

import re

times = []
memory = []

for stats in glob.glob('*.stats'):
    with open(stats) as f:
        lines = f.readlines()
        string = ''.join(lines)

        mem = re.findall(r'\d+K memory used', string)
        mem_string = ''.join(mem)
        mem_d = re.findall(r'\d+', mem_string)
        try:
            memory.append(int(mem_d[0])/1000000)
        except:
            print(stats)

        time = re.findall(r'\d+:\d+:\d+ walltime', string)
        time_string = ' '.join(time)
        time_d = re.findall(r'\d+:\d+:\d+', time_string)
        try:
            timestamp = datetime.strptime(time_d[1], '%H:%M:%S')
            mins = (timestamp - datetime(1900, 1, 1, 0, 0, 0)).seconds / 60
            times.append(mins)
        except:
            print(stats)


plt.figure()
sb.set_context("talk")
plt.title("Memory taken to run queries")
sb.displot(x=memory, bins=30)
sb.despine()
plt.ylabel("Frequency")
plt.xlabel("Memory used (GB)")
plt.xticks(rotation=30, ha='right')
plt.savefig("Test_mem.png",bbox_inches='tight')

plt.figure()
sb.set_context("talk")
plt.title("Time taken to run queries")
sb.displot(x=times, bins=30)
sb.despine()
plt.ylabel("Frequency")
plt.xlabel("Time taken (mins)")
plt.xticks(rotation=30, ha='right')
plt.savefig("Test_times.png",bbox_inches='tight')
