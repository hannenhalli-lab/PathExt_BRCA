import pandas as pd
from collections import defaultdict
import statsmodels.stats.multitest as bh
import random
random.seed(1)
import networkx as nx
import numpy as np
import sys
import math

import microarray_functions as mic_fun
import network_functions as net_fun
import percentile_functions as perc_fun

if len(sys.argv) != 7:
	print("argv[1] = node weights fname")
	print("argv[2] = colname of condition to study")
	print("argv[3] = unweighted network file")
	print("argv[4] = path length threshold")
	print("argv[5] = file with paths, z-scores and p-values (before FDR)")
	print("argv[6] = output file name")
	sys.exit(1)

# Set user inputs
node_weight_fname = sys.argv[1]
condition = sys.argv[2] # This is the condition we want to study
unweighted_nw_fname = sys.argv[3]
path_length_thresh = int(sys.argv[4])
zscore_pvals_fname = sys.argv[5]
out_fname = sys.argv[6]

# Set other inputs
percentiles = list(np.arange(0.001, 0.01, 0.001))
percentiles.extend(list(np.arange(0.01, 0.1, 0.01)))
percentiles.extend(list(np.arange(0.1, 0.6, 0.1)))
qscores = [0.05, 0.01, 0.005, 0.001]
default_alpha = 0.05

# Read z-scores and corresponding p-values
zscore_pvals = pd.read_table(zscore_pvals_fname)
zscore_pvals.columns = ('ij', 'normaltest_before_boxcox', 'normaltest_after_boxcox', 'zscore_after_boxcox', 'two_sided_pval_for_zscore')
zscore_pvals = zscore_pvals.set_index(zscore_pvals.columns.values[0])
zscore_pvals = zscore_pvals.dropna()
print("Done reading z-scores and corresponding p-vals. Shape is ", zscore_pvals.shape)

# Read node weights
# Column 0 -> gene labels. Make this the index
# Columns 1 through m -> node weights for the various conditions
node_weight = pd.read_csv(node_weight_fname, sep = '\t', index_col = 0)
node_weight = node_weight.groupby(by=node_weight.index).first() # if there are multiple rows with same index, keep first row
node_weight = mic_fun.restructure_node_weight(node_weight, condition) # column 0 -> condition of interest, 2 to m -> other conditions
print("Read node weight matrix with ", node_weight.shape[0], " rows and ", node_weight.shape[1], " columns")
#print(node_weight.head())

# Read unweighted network
G_unweighted = nx.read_edgelist(unweighted_nw_fname, delimiter = "\t", nodetype = str, create_using = nx.DiGraph())
print("Read network with ", len(G_unweighted.nodes()), " nodes and ", len(G_unweighted.edges()), " edges")

# Map node weights onto unweighted network
node_weight_relevant = mic_fun.get_relevant_node_weight(node_weight)
G_base = net_fun.get_weighted_network(node_weight_relevant, G_unweighted)
print("Got base network with ", len(G_base.nodes()), " nodes and ", len(G_base.edges()), " edges")

# Get all-pairs-shortest-paths
Pij = net_fun.get_all_sp_paths_costs(G_base)
print("Got shortest path costs for ", Pij.shape[0], " node-pairs")

# For each percentile and each q-score threshold, figure out the size of the top-net
# For a percentile, get the paths within that. Get the corresponding z-scores and p-vals
# Then within that, apply FDR
# Then try different q-score thresholds for the top-net
with open(out_fname, 'w') as f:
	f.write( '\t'.join(['percentile', 'qscore', 'topnet_nodes', 'topnet_edges']) + '\n' )
	for percentile in percentiles:
		print("")
		print("############### Percentile ", percentile, " ###############")

		# Get paths within this percentile cutoff (accounting for path length threshold)
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

		for qscore in qscores:
			print("############### q-score ", qscore, " ###############")
			significant_paths = zscore_pvals_temp.loc[zscore_pvals_temp['bh_qscore'] <= qscore]
			print("Got ", significant_paths.shape[0]," paths within ", percentile, " percentile and q-score <= ", qscore)
			print(significant_paths.head())

			# Extract top-net using the significant paths
			topnet_edges = set()
			topnet_nodes = set()
			for path_str in significant_paths.index:
				path = path_str.split('#')
				topnet_nodes = topnet_nodes.union(set(path))
				topnet_edges = topnet_edges.union(net_fun.get_edges_in_path(path, G_base))
			print("Got ", len(topnet_nodes), " nodes and ", len(topnet_edges), " edges in top-net")

			f.write( '\t'.join([str(percentile), str(qscore), str(len(topnet_nodes)), str(len(topnet_edges))]) + '\n' )
