#!/bin/bash
kin_mtx_block_path="/Data/Kinship_matrix/freeze.8.JHS_mtx/JHS_kinship_matrix_from_dense_TOPMed_freeze_8.rds"
kin_mtrx_path="/Data/Kinship_matrix/freeze.8.JHS_mtx_3X/JHS_kinship_matrix_from_dense_TOPMed_freeze_8_3X.rds"

output_path="/Results/performance_measure/3X/"
bin_size=0.1
rho_e=0 
n_sim=10000
num_blocks=3 #or 10X

#paths for results of performance measures
res_df_mvnfast_by_block_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_mvnfast_by_blocks.RData"
res_df_mvnfast_over_all_blocks_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_mvnfast_over_all_blocks.RData"

res_df_calcGC_KMZ_by_block_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcGC_KMZ_by_block_path.RData"
res_df_calcGC_KMA_by_block_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcGC_KMA_by_block_path.RData"
res_df_calcHE_1_by_block_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcHE_1_by_block_path.RData"
res_df_calcHE_2_by_block_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcHE_2_by_block_path.RData"

res_df_calcGC_KMZ_final_over_all_blocks_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcGC_KMZ_over_all_blocks_path.RData"
res_df_calcGC_KMA_final_over_all_blocks_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcGC_KMA_over_all_blocks_path.RData"
res_df_calcHE_1_final_over_all_blocks_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcHE_1_over_all_blocks_path.RData"
res_df_calcHE_2_final_over_all_blocks_path="/Code/performance_measure/Outputs/3X/${job_id}_res_df_calcHE_2_over_all_blocks_path.RData"


for job_id in {1..10} 
do
   3_parametric_bootstrap.R $job_id $kin_mtx_block_path \
      $kin_mtrx_path $output_path $bin_size $rho_e $n_sim $num_blocks \
      $res_df_mvnfast_by_block_path $res_df_mvnfast_over_all_blocks_path \
      $res_df_calcGC_KMZ_by_block_path $res_df_calcGC_KMA_by_block_path \
      $res_df_calcHE_1_by_block_path $res_df_calcHE_2_by_block_path \
      $res_df_calcGC_KMZ_final_over_all_blocks_path $res_df_calcGC_KMA_final_over_all_blocks_path \
      $res_df_calcHE_1_final_over_all_blocks_path $res_df_calcHE_2_final_over_all_blocks_path \

done