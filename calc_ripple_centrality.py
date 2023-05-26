import sys
import networkx as nx

# Print first n key, value pairs
def print_first_n(some_dict, n):
	num = 0
	for key, value in some_dict.items():
		if num > n:
			return
		print(key, value)
		num += 1

def outward_reachability(graph):
	shortest_path_costs = dict(nx.all_pairs_dijkstra_path_length(graph))
	reachable = {} # key -> node, value -> no. of nodes reachable from that node
	for n in graph.nodes():
		reachable[n] = len(shortest_path_costs[n].keys())
	return reachable

def ripple_centrality(closeness_centr, outward_reach, num_nodes):
	ripple_centrality = {}
	for node, closeness in closeness_centr.items():
		ripple_centrality[node] = closeness * (outward_reach[node]/num_nodes)
	return ripple_centrality


if len(sys.argv) != 3:
	print("argv[1] = weighted network file")
	print("argv[2] = output file")
	sys.exit(1)

# read network file
G = nx.read_edgelist(sys.argv[1], nodetype=str, data=(('weight',float),), create_using = nx.DiGraph())
print("Read network")
print("No. of edges = ", len(G.edges()))
print("No. of nodes = ", len(G.nodes()))

############### Calculate centralities ##################
# in degree centrality
#in_degree_centr = nx.in_degree_centrality(G)
#print("Calculated in-degree centrality")

# out degree centrality
#out_degree_centr = nx.out_degree_centrality(G)
#print("Calculated out-degree centrality")

# closeness centrality
closeness_centr = nx.closeness_centrality(G, distance = 'weight')
print("Calculated closeness centrality")

# betweenness centrality
#betweenness_centr = nx.betweenness_centrality(G, weight='weight')
#print("Calculated betweenness centrality")

# eigenvector centrality
#eigenvector_centr = nx.eigenvector_centrality_numpy(G, weight='weight')
#print("Calculated eigenvector centrality")

# outward reachability
outward_reach = outward_reachability(G)
print("Calculated outward reachability")

# ripple centrality
ripple_centr = ripple_centrality(closeness_centr, outward_reach, len(G.nodes()))
print("Calculated ripple centrality")

############### Print centralities ##################
with open(sys.argv[2], "w") as f:
	f.write("node" + "\t")
	#f.write("indegree" + "\t")
	#f.write("indegree_centrality" + "\t")
	#f.write("outdegree" + "\t")
	#f.write("outdegree_centrality" + "\t")
	#f.write("betweenness_centrality" + "\t")
	#f.write("eigenvector_centrality" + "\t")
	f.write("closeness_centrality" + "\t")
	f.write("outward_reachability" + "\t")
	f.write("ripple_centrality" + "\n")
	for node in G.nodes():
		f.write(node + "\t")
		#f.write(str(G.in_degree(node)) + "\t")
		#f.write(str(in_degree_centr[node]) + "\t")
		#f.write(str(G.out_degree(node)) + "\t")
		#f.write(str(out_degree_centr[node]) + "\t")
		#f.write(str(betweenness_centr[node]) + "\t")
		#f.write(str(eigenvector_centr[node]) + "\t")
		f.write(str(closeness_centr[node]) + "\t")
		f.write(str(outward_reach[node]) + "\t")
		f.write(str(round(ripple_centr[node],4)) + "\n")
