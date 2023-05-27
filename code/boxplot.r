### Python code ####
#### Boxplot GC distribution #####

import sys
import pandas as pd
import numpy as np
import seaborn as sns
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

#### Read input file ###
df = pd.read_csv('boxplot_input.csv')
order = ['PathExt','DEGs']

sns.set_style('whitegrid')
ax= sns.boxplot(data=df, order=order, notch=True)
ax= sns.stripplot(data=df, order=order)
ax.set(ylim=(0,1))

add_stat_annotation(ax, data=df, order=order, box_pairs=[("PathExt", "DEGs")], test='Mann-Whitney', text_format='full', loc='outside', verbose=2)

#ax.set_ylabel("Depmap Probability Score",fontsize=15,weight='bold')
ax.set_xlabel("Basal",fontsize=15,weight='bold')
ax.figure.set_size_inches(5,5)

plt.savefig('Boxplot_out', dpi=600i)

