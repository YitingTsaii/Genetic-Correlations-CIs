library(plyr)
library(ggplot2)
library(EnvStats)


source("parametric_bootstrap_helper_functions.R")
source("combine_parametric_bootstrap_simulations_function.R")
source("compute_CI_functions.R")
source("stats_analysis_functions.R")



args <- commandArgs(trailingOnly=TRUE)

simulation_database_path <- args[1]
simulation_data_useful_format <- args[2]

bin_size <- as.numeric(args[3])
alpha <- as.numeric(args[4])

cover_prob_all_path_pmf <- args[5]
CI_width_all_path_pmf <- args[6]

cover_prob_settingA_pmf_path <- args[7]
CI_width_settingA_pmf_path <- args[8]

# (b) compute CI
# input: h1_sq_hat, h2_sq_hat, rho_k_hat, bin_size, alpha, path_to_combined_results
load(simulation_data_useful_format) 
pmf_adj <- parametric_bootstrap_combined_results$pmf_adj
pmf_zerodiag <- parametric_bootstrap_combined_results$pmf_zerodiag


# get cover_prob_all, CI_width_all
# for pmf 
# run for all different kinship setting (i.e. setting A (3X), setting B (10X))
dfs_pmf <- get_summary_df_pmf(simulation_database_path, pmf_adj, pmf_zerodiag) 
cover_prob_all_pmf <- dfs_pmf$cover_prob_all
CI_width_all_pmf <- dfs_pmf$CI_width_all

save(cover_prob_all_pmf,file=cover_prob_all_path_pmf) 
save(CI_width_all_pmf,file=CI_width_all_path_pmf) 


# create and save the data frames for plotting 
a <- apply(cover_prob_all_pmf[,2:101],1,mean)
cover_prob_settingA_pmf <- data.frame(true_rho = cover_prob_all_pmf$true_rho_k, cover_prob = a)  
save(cover_prob_settingA_pmf, file=cover_prob_settingA_pmf_path)

b <- apply(CI_width_all_pmf[,2:101],1,mean)
CI_width_settingA_pmf <- data.frame(true_rho = CI_width_all_pmf$true_rho_k, CI_width = b)  
save(CI_width_settingA_pmf, file=CI_width_settingA_pmf_path)
