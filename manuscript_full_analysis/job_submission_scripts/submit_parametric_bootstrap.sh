#!/bin/bash
#SBATCH --job-name="param_boot"    # Job name
#SBATCH --ntasks=1             # Run a single task
#SBATCH --time=96:00:00        # Time limit hrs:min:sec
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --output="/Code/logs_and_errors_matrix_10X/param_boot_%A_%a.log"
#SBATCH --error="/Code/logs_and_errors_matrix_10X/param_boot_%A_%a.err"
#SBATCH --array=[1-2000] 


#Note: repeat for 3X
job_id=$SLURM_ARRAY_TASK_ID

kin_mtx_block_path="/Data/Kinship_matrix/freeze.8.JHS_mtx/JHS_kinship_matrix_from_dense_TOPMed_freeze_8.rds"
kin_mtrx_path="/Data/Kinship_matrix/freeze.8.JHS_mtx_10X/JHS_kinship_matrix_from_dense_TOPMed_freeze_8_10X.rds"

output_path="/Results/bootstrap_out_matrix_10X/"

bin_size=0.1
rho_e=0
n_sim=10000
num_blocks=10 

Rscript 2_parametric_bootstrap.R $job_id $kin_mtx_block_path $kin_mtrx_path $output_path $bin_size $rho_e $n_sim $num_blocks



