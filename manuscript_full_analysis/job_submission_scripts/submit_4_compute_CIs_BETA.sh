#!/bin/bash
#SBATCH --job-name="CI_BETA"    
#SBATCH --ntasks=1             
#SBATCH --time=96:00:00       
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --output="/Code/logs_and_errors/compute_CI/CI_BETA_%A_%a.log"
#SBATCH --error="/Code/logs_and_errors/compute_CI/CI_BETA_%A_%a.err"


#Note: repeat for 10X

simulation_database_path="/Results/data_simulations/3X"

simulation_data_useful_format="/Results/combined_data_3x/3X_matrix_parametric_bootsrtap.RData"


bin_size=0.1

alpha=0.05

cover_prob_all_path_beta="/Results/CI_results/3X/cover_prob_all/cover_prob_all_setting_3X_beta.RData"

CI_width_all_path_beta="/Results/CI_results/3X/CI_width_all_path/CI_width_all_setting_3X_beta.RData"

cover_prob_settingA_beta_path="/Results/CI_results/3X/files_for_plotting/cover_prob_all_setting_3X_beta.RData"

CI_width_settingA_beta_path="/Results/CI_results/3X/files_for_plotting/CI_width_all_setting_3X_beta.RData"

Rscript 4_compute_CI_from_parametric_bootstrap_BETA.R $simulation_database_path $simulation_data_useful_format $bin_size $alpha $cover_prob_all_path_beta $CI_width_all_path_beta $cover_prob_settingA_beta_path $CI_width_settingA_beta_path


