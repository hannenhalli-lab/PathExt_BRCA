#!/bin/sh

trim_last_slash () {
        wdir=$1
        if [ "${wdir: -1}" == "/" ]
        then
                wdir="${wdir%?}"
        fi
        echo "${wdir}"
}

if [ "$#" -ne 5 ]
then
	echo "argv[1] = node weight file"
	echo "argv[2] = colname of condition to study"
	echo "argv[3] = unweighted network file"
	#echo "argv[5] = percentile threshold"
	#echo "argv[6] = path length threshold"
	#echo "argv[7] = q-score cutoff"
	#echo "argv[8] = number of randomizations"
	echo "argv[4] = output directory"
	echo "argv[5] = file name for base network (we'll put it in the output directory)"
	#echo "argv[11] = file name for TopNet (we'll put it in the output directory)"
        exit
fi

# Set inputs
data_fname=${1}
condition_of_interest=${2}
unweighted_nw_fname=${3}
percentile=0.1
path_length_thresh=2
#qscore_thresh=${7}
num_trials=1000
out_dir=$(trim_last_slash ${4})
base_nw_fname=${5}
#topnet_fname=${11}

# Create a temporary directory in the output directory
# We will put all the Pij files here
# Once the z-score has been calculated we will remove the Pij files as well as the temp directory
mkdir ${out_dir}/temp

# Map the node weights onto the unweighted network to get an base network.
# Compute top 'percentile' shortest paths in this (actual) network. Write these paths and costs.
# Randomize the node weights 'num_trials' times, resulting in 'num_trials' randomized networks.
# Compute the cost of the same paths as in the top 'percentile' shortest paths in the actual network.
# Write these paths and costs into output files.
python node_weight_matrix_colname_Pijs.py ${data_fname} ${condition_of_interest} ${unweighted_nw_fname} ${percentile} ${path_length_thresh} ${num_trials} ${out_dir}/${base_nw_fname} ${out_dir}/temp/Pij

# Calculate z-score and corresponding p-value for each path
python fdr_rand_pijs_boxcox.py ${out_dir}/temp ${out_dir}/Pij_zscores.txt

# Delete all temporary files
rm -rf ${out_dir}/temp

# Run try different thresholds
python try_different_thresholds_node_weight_matrix.py ${data_fname} ${condition_of_interest} ${unweighted_nw_fname} ${path_length_thresh} ${out_dir}/Pij_zscores.txt ${out_dir}/thresh_TopNet_sizes.txt

# Carry out FDR (benjamini-hochberg)
#python benjamini_hochberg_boxcox.py ${out_dir}/Pij_zscores.txt ${qscore_thresh} ${out_dir}/Pij_zscores_fdr.txt

# Extract top-net based on paths whose q-score <= qscore_thresh
#python extract_fdr_network.py ${out_dir}/${base_nw_fname} ${out_dir}/Pij_zscores_fdr.txt ${qscore_thresh} ${out_dir}/${topnet_fname}
