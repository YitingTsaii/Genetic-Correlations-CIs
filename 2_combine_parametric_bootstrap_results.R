library(plyr)
library(ggplot2)
library(EnvStats)


source("parametric_bootstrap_helper_functions.R")
source("combine_parametric_bootstrap_simulations_function.R")
source("compute_CI_functions.R")
source("stats_analysis_functions.R")


parametric_bootstrap_database_path <- "/Results/parametric_bootstrap/simulations_out/"
bin_size <- 0.1
out_path <- "/Results/parametric_bootstrap/simulations_out/combined/combined_parametric_bootsrtap_results.RData"



# re-format results and save 
parametric_bootstrap_combined_results <- combine_simulation_data(parametric_bootstrap_database_path, bin_size)
save(parametric_bootstrap_combined_results,file = out_path)
