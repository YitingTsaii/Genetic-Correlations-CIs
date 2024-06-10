library(plyr)
library(ggplot2)
library(EnvStats)


source("parametric_bootstrap_helper_functions.R")
source("combine_parametric_bootstrap_simulations_function.R")
source("compute_CI_functions.R")
source("stats_analysis_functions.R")


parametric_bootstrap_database_path <- "/Results/parametric_bootstrap/10X_matrix/simulations_out_matrix_10X_1_2000/simulations_out_matrix_10X_1_to_2000/simulations_out_matrix_10X/"
bin_size <- 0.1
out_path <- "/Results/parametric_bootstrap/10X_matrix/simulations_out_matrix_10X_1_2000/simulations_out_matrix_10X_1_to_2000/combined_data_10X_complete/10X_matrix_parametric_bootsrtap.RData"



# re-format results and save 
parametric_bootstrap_combined_results <- combine_simulation_data(parametric_bootstrap_database_path, bin_size)
save(parametric_bootstrap_combined_results,file = out_path)
