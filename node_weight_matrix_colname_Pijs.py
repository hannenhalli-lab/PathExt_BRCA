# Code to create a network given inputs PPI and node weight matrix,
# and track top x percentile shortest paths in the actual data
# We also create a background of randomized data by permuting the input node weights
# For each resulting randomized network, we compute the cost of the paths that are in the
# top x percentile in the actual data

# NOTE: This code assumes the node weight data has only 1 column per condition
#	Column 0 -> Gene symbols, 
# 	Columns 1 through m -> node weight in the conditions under study
#	Header for the condition of interest is given as user input

import pandas as pd
import random
import networkx as nx
import numpy as np
import sys
import math

import microarray_functions as mic_fun
import network_functions as net_fun
import percentile_functions as perc_fun

# node weights are given as a pandas dataframe, indexed by gene
# column 0 -> node weight in condition of interest
def combine_data_get_sp_paths_costs_basenw(G_ref, node_weight, base_nw_fname, paths):
	node_weight_relevant = mic_fun.get_relevant_node_weight(node_weight)
	print("Relevant node weight matrix has shape ", node_weight_relevant.shape)

	# Map node weights onto reference network
	G_base = net_fun.get_weighted_network(node_weight_relevant, G_ref)
	print("Got base network with ", len(G_base.nodes()), " nodes and ", len(G_base.edges()), " edges")

	# Drop node_weight values for genes which don't map to base network
	genes_to_drop = set(node_weight.index) - set(G_base.nodes())
	node_weight = node_weight.drop(genes_to_drop)

	if len(paths) == 0: # Actual data
		# Write the base network to file
		nx.write_weighted_edgelist(G_base, base_nw_fname, delimiter = '\t')
		# Get all-pairs-shortest-path costs
		# Return value is a pandas dataframe, indexed by the string src#dest
		Pij = net_fun.get_all_sp_paths_costs(G_base)
		print("Got shortest path costs for ", Pij.shape[0], " node-pairs")
	elif len(paths) > 0: # Randomized data
		# Get cost of the same paths that are the shortest in the actual dataset
		Pij = net_fun.get_costs_of_given_paths(G_base, paths)
		print("Got cost of paths which are shortest in the actual data")

	return Pij, node_weight


if len(sys.argv) != 9:
	print("argv[1] = node weight file")
	print("argv[2] = colname of condition to study")
	print("argv[3] = reference network (weighted edgelist)")
	print("argv[4] = percentile threshold")
	print("argv[5] = path length threshold")
	print("argv[6] = number of randomizations")
	print("argv[7] = output file for base network")
	print("argv[8] = output file prefix for Pij (including directory)")
	sys.exit(1)

#SMALL_VAL = pow(10, -6)

# Set inputs
node_weight_fname = sys.argv[1]
condition = sys.argv[2] # This is the condition we want to study
ref_nw_fname = sys.argv[3]
percentile = float(sys.argv[4]) # We'll only keep paths whose cost < this threshold
path_length_thresh = int(sys.argv[5]) # We'll only keep paths with length >= this threshold
num_trials = int(sys.argv[6])
base_nw_fname = sys.argv[7]
output_fname_prefix = sys.argv[8]

# Read base node_weight
# Column 0 -> gene labels. Make this the index
# Columns 1 through m -> node weights in various conditions
node_weight = pd.read_csv(node_weight_fname, sep = '\t', index_col = 0)
node_weight = node_weight.groupby(by=node_weight.index).first() # if there are multiple rows with same index, keep first row
#node_weight = node_weight + SMALL_VAL

node_weight = mic_fun.restructure_node_weight(node_weight, condition) # column 0 -> condition of interest, 2 to m -> other conditions
print("Read node weight matrix with ", node_weight.shape[0], " rows and ", node_weight.shape[1], " columns")
#print(node_weight.head())

# Read reference network
G_ref = nx.read_weighted_edgelist(ref_nw_fname, delimiter = "\t", nodetype = str, create_using = nx.DiGraph())
print("Read network with ", len(G_ref.nodes()), " nodes and ", len(G_ref.edges()), " edges")

# Compute Pij in actual data
Pij, node_weight = combine_data_get_sp_paths_costs_basenw(G_ref, node_weight, base_nw_fname, []) # Give an empty list as the paths argument
print("After dropping genes which don't map to our network, got node weights for ", node_weight.shape[0], " genes and ", node_weight.shape[1], " conditions")
print(node_weight.head())

Pij = perc_fun.get_Pij_percentile_norm_cost(Pij, percentile, path_length_thresh)
print("After taking percentile cutoff, Pij has ", Pij.shape[0], " rows and ", Pij.shape[1], " columns")
print(Pij.head())
Pij.to_csv(output_fname_prefix+"_actual.txt", sep = "\t")

# Randomized data
for i in range(num_trials):
	print("######################## Trial ", i, " ########################")

	randomized_node_weight = mic_fun.get_randomized_node_weight_matrix(node_weight)
	print("After shuffling, got node weight values for ", randomized_node_weight.shape[0], " genes and ", randomized_node_weight.shape[1], " conditions")
	print(randomized_node_weight.head())

	# Give shortest paths in actual dataset as the paths argument
	Pij, randomized_node_weight = combine_data_get_sp_paths_costs_basenw(G_ref, randomized_node_weight, base_nw_fname, list(Pij.index))

	Pij.to_csv(output_fname_prefix+"_"+str(i)+".txt", sep = "\t")
