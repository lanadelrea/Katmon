#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

from os import rename

# Input
position_base_counts_csv = sys.argv[1]
position_base_counts = pd.read_csv(position_base_counts_csv, sep=',')

bammix_plot = sys.argv[2]

sample = sys.argv[3]

# Create a plot from the positions and proportions
df_bam = position_base_counts[["Position", "A_proportion", "C_proportion", "G_proportion","T_proportion"]].copy()
df_bam.rename(columns = {'A_proportion':'A', 'C_proportion':'C', 'G_proportion':'G', 'T_proportion': 'T'}, inplace=True)
df_bam.set_index("Position", inplace=True)
#df_bam = df_bam[~(df_bam == 0).all(axis=1)]

df_bam.plot(kind='bar', stacked=True, color=[ "#186F65", "#B2533E", "#7C81AD", "#EE9322"], figsize=(9, 5))

# Labels for x and y axis
plt.xlabel('Position')
plt.ylabel('Proportion')

# Other plot stuff
plt.title(sample)

plt.xticks(rotation=90, fontsize=5)
plt.ylim([0,1])
plt.legend(bbox_to_anchor=[1,1], ncol=1)
plt.tight_layout()
plt.savefig(f"{sample}_bammix_plot", dpi=200) # Name plot file according to sample