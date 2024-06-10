library(plyr)
library(ggplot2)
library(EnvStats)



source("parametric_bootstrap_helper_functions.R")
source("combine_parametric_bootstrap_simulations_function.R")
source("compute_CI_functions.R")
source("stats_analysis_functions_for_FISHER.R")


args <- commandArgs(trailingOnly=TRUE)

simulation_database_path <- args[1]
simulation_data_useful_format <- args[2]

bin_size <- as.numeric(args[3])
alpha <- as.numeric(args[4])

cover_prob_all_path_fisher <- args[5]
CI_width_all_path_fisher <- args[6]

cover_prob_settingA_fisher_path <- args[7]
CI_width_settingA_fisher_path <- args[8]

kinMat_zero_diag_path <- args[9]


kinMat.zero.diag <- readRDS(kinMat_zero_diag_path)

kinMat.adj_n <- sqrt(sum(diag(kinMat.zero.diag %*% t(kinMat.zero.diag))))


# get cover_prob_all, CI_width_all
# for fisher
# run for all different kinship setting (i.e. setting A (3X), setting B (10X))

dfs_fisher <- get_summary_df_fisher(simulation_database_path, kinMat.adj_n)
cover_prob_all_fisher <- dfs_fisher$cover_prob_all
CI_width_all_fisher <- dfs_fisher$CI_width_all

save(cover_prob_all_fisher,file=cover_prob_all_path_fisher) 
save(CI_width_all_fisher,file=CI_width_all_path_fisher) 


# create and save the data frames for plotting
a <- apply(cover_prob_all_fisher[,2:101],1,mean)
cover_prob_settingA_fisher <- data.frame(true_rho = cover_prob_all_fisher$true_rho_k, cover_prob = a) # change the variable name according to the setting (A, B, C) and the method (pmf, beta, fisher)
save(cover_prob_settingA_fisher, file=cover_prob_settingA_fisher_path)


b <- apply(CI_width_all_fisher[,2:101],1,mean)
CI_width_settingA_fisher <- data.frame(true_rho = CI_width_all_fisher$true_rho_k, CI_width = b) # change the variable name according to the setting (A, B, C) and the method (pmf, beta, fisher) 
save(CI_width_settingA_fisher, file=CI_width_settingA_fisher_path)

