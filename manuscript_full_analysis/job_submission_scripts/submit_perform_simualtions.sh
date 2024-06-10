#!/bin/bash
#SBATCH --job-name="ps_3X"    # Job name
#SBATCH --ntasks=1             # Run a single task
#SBATCH --time=96:00:00        # Time limit hrs:min:sec
#SBATCH --cpus-per-task=3
#SBATCH --exclusive
#SBATCH --output="/Code/logs_and_errors/sim_data_%A_%a.log"
#SBATCH --error="/Code/logs_and_errors/sim_data_%A_%a.err"
#SBATCH --array=[1-2000] 

#note: repeat for 10X

job_id=$SLURM_ARRAY_TASK_ID

kin_mtx_block_path="/Data/Kinship_matrix/freeze.8.JHS_mtx/JHS_kinship_matrix_from_dense_TOPMed_freeze_8.rds"
kin_mtrx_path="/Data/Kinship_matrix/freeze.8.JHS_mtx_3X/JHS_kinship_matrix_from_dense_TOPMed_freeze_8_3X.rds"

output_path="/Results/data_sims/simulations_out_matrix_3X/"
bin_size=0.1
rho_e=0 
n_sim=1000 
num_blocks=3 #10X


Rscript 3_perform_simulations.R $job_id $kin_mtx_block_path $kin_mtrx_path $output_path $bin_size $rho_e $n_sim $num_blocks