import pandas as pd
from collections import defaultdict
import statsmodels.stats.multitest as bh
import random
import networkx as nx
import numpy as np
import sys
import math

import microarray_functions as mic_fun
import network_functions as net_fun
import percentile_functions as perc_fun

if len(sys.argv) != 10:
	print("argv[1] = microarray data fname")
	print("argv[2] = name of perturbation sample to study")
	print("argv[3] = name of control sample")
	print("argv[4] = unweighted network file")
	print("argv[5] = path length threshold")
	print("argv[6] = percentile threshold")
	print("argv[7] = q-score threshold")
	print("argv[8] = file with paths, z-scores and p-values (before FDR)")
	print("argv[9] = output file name for topnet")
	sys.exit(1)

# Set user inputs
data_fname = sys.argv[1]
perturbation_sample = sys.argv[2] # This is the perturbation we want to study
control_sample = sys.argv[3]
unweighted_nw_fname = sys.argv[4]
path_length_thresh = int(sys.argv[5])
percentile = float(sys.argv[6])
qscore = float(sys.argv[7])
zscore_pvals_fname = sys.argv[8]
out_fname = sys.argv[9]

# Set other inputs
#percentiles = [0.1, 0.2, 0.3, 0.4, 0.5]
#qscores = [0.05, 0.01, 0.005, 0.001]
default_alpha = 0.05

# Read z-scores and corresponding p-values
zscore_pvals = pd.read_table(zscore_pvals_fname)
zscore_pvals.columns = ('ij', 'normaltest_before_boxcox', 'normaltest_after_boxcox', 'zscore_after_boxcox', 'two_sided_pval_for_zscore')
zscore_pvals = zscore_pvals.set_index(zscore_pvals.columns.values[0])
zscore_pvals = zscore_pvals.dropna()
print("Done reading z-scores and corresponding p-vals. Shape is ", zscore_pvals.shape)

# Read microarray data
# Column 0 -> gene labels. Make this the index
# Columns 1 through m -> a control, and various perturbed conditions
SI = pd.read_table(data_fname, na_values=' ')
SI = SI.set_index(SI.columns[0])
SI = mic_fun.restructure_SI(SI, perturbation_sample, control_sample) # column 0 -> perturbation to study, 
                                                                     # column 1 -> control, 
                                                                     # columns 2 to m -> other perturbations
print("Read microarray data with ", SI.shape[0], " rows and ", SI.shape[1], " columns")

# Read unweighted network
G_unweighted = nx.read_edgelist(unweighted_nw_fname, delimiter = "\t", nodetype = str, create_using = nx.DiGraph())
print("Read network with ", len(G_unweighted.nodes()), " nodes and ", len(G_unweighted.edges()), " edges")

# Map gene expression values onto unweighted network
SI_relevant = mic_fun.get_relevant_SI(SI)
G_response = net_fun.get_repressed_response_network(SI_relevant, G_unweighted)
print("Got response network with ", len(G_response.nodes()), " nodes and ", len(G_response.edges()), " edges")

# Get all-pairs-shortest-paths
Pij = net_fun.get_all_sp_paths_costs(G_response)
print("Got shortest path costs for ", Pij.shape[0], " node-pairs")

####### For the given percentile and q-score threshold, figure out the size of the top-net
####### Get the paths within that percentile threshold. Get the corresponding z-scores and p-vals
####### Then within that, apply FDR. Then pick only paths within the q-score threshold
####### Get paths within this percentile cutoff (accounting for path length threshold)
Pij_temp = perc_fun.get_Pij_percentile_norm_cost(Pij, percentile, path_length_thresh)
print("After taking percentile cutoff, Pij has ", Pij_temp.shape[0], " rows and ", Pij_temp.shape[1], " columns")

# Get z-scores and corresponding p-values for these paths
zscore_pvals_temp = zscore_pvals.join(Pij_temp, how = 'inner')
print("After taking inner join, shape of zscore_pvals_temp is ", zscore_pvals_temp.shape)
print(zscore_pvals_temp.head())

# Apply FDR only for this set of paths
bh_output = bh.multipletests(zscore_pvals_temp['two_sided_pval_for_zscore'], alpha = default_alpha, method = 'fdr_bh')
reject = bh_output[0]
bh_pval = bh_output[1]
zscore_pvals_temp['bh_qscore'] = bh_pval
print("Done with BH test")
print(zscore_pvals_temp.head())

significant_paths = zscore_pvals_temp.loc[zscore_pvals_temp['bh_qscore'] <= qscore]
print("Got ", significant_paths.shape[0]," paths within ", percentile, " percentile and q-score <= ", qscore)
print(significant_paths.head())

# Extract top-net using the significant paths
topnet_edges = set()
topnet_nodes = set()
for path_str in significant_paths.index:
	path = path_str.split('#')
	topnet_nodes = topnet_nodes.union(set(path))
	topnet_edges = topnet_edges.union(net_fun.get_edges_in_path(path, G_response))
print("Got ", len(topnet_nodes), " nodes and ", len(topnet_edges), " edges in top-net")

# Write top-net
with open(out_fname, 'w') as f:
	for edge in topnet_edges:
		f.write('\t'.join([edge[0], edge[1], str(edge[2])]) + '\n')
