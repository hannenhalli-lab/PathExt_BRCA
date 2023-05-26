#!/bin/bash

module load python
mkdir test_data/results/temp
python node_weight_matrix_colname_Pijs.py test_data/input_data Mean test_data/human_PPIN.txt 0.1 2 1000 test_data/results/Activated_response test_data/results/temp/Pij
python /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/fdr_rand_pijs_boxcox.py /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/results/temp /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/results/Pij_zscores.txt
rm -rf /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/results/temp
python /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/try_different_thresholds_node_weight_matrix.py /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/input_data Mean /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/human_PPIN.txt 2 /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/results/Pij_zscores.txt /data/agrawalp4/covid_project/ICU_nonICU/icu_mean_active/test_data/results/thresh_TopNet_sizes.txt

