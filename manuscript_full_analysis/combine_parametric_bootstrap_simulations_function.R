######### functions to combine simulation data into a useful format #############
combine_simulation_data <- function(simulation_database_path, bin_size){
  # get sum_h1h2 (red part) -> 20 of them
  
  sum_h1h2_adj <- list()
  sum_h1h2_zerodiag <- list()
  i <- 1
  start <- -(1+(1-bin_size))/2
  
  end <- (1+(1-bin_size))/2
 
  for(true_rho_k in seq(start, end, by=bin_size)){ # loop over true_rho_k (20 in total)
    true_rho_k <- round(true_rho_k, 3)
    file_names <- list.files(path=simulation_database_path, pattern = paste0("rhok-", true_rho_k)) # file_names should be 10*10 in length
    sum_h1h2_each_adj <- array(0, dim = c(10,10,20))
    sum_h1h2_each_zerodiag <- array(0, dim = c(10,10,20))
    for(file_name in file_names){ # 10*10 (loop over h1, h2)
      load(paste0(simulation_database_path, file_name)) #rres
      bins_adj_prob <- rres$bins_3D_adj/sum(rres$bins_3D_adj)  # check: sum(bins_adj_prob) = 1
      bins_zerodiag_prob <- rres$bins_3D_zerodiag/sum(rres$bins_3D_zerodiag)  # check: sum(bins_zerodiag_prob) = 1
      sum_h1h2_each_adj <- sum_h1h2_each_adj + bins_adj_prob
      sum_h1h2_each_zerodiag <- sum_h1h2_each_zerodiag + bins_zerodiag_prob
    }
    sum_h1h2_adj[[i]] <- sum_h1h2_each_adj
    sum_h1h2_zerodiag[[i]] <- sum_h1h2_each_zerodiag
    i <- i+1
    print(i)
  }
  
  # get sum_h1h2rhok (blue part) -> only 1 result
  sum_h1h2rhok_adj <- array(0, dim = c(10,10,20))
  sum_h1h2rhok_zerodiag <- array(0, dim = c(10,10,20))
  for(i in 1:20){
    sum_h1h2rhok_adj <- sum_h1h2rhok_adj + sum_h1h2_adj[[i]]
    sum_h1h2rhok_zerodiag <- sum_h1h2rhok_zerodiag + sum_h1h2_zerodiag[[i]]
  }
  
  # get the new pmf
  pmf_adj <- list()
  pmf_zerodiag <- list()
  for(i in 1:20){ # i is the index of true rho_k
    pmf_adj[[i]] <- sum_h1h2_adj[[i]] / sum_h1h2rhok_adj # element-wise division # a 10*10*20 matrix divided by a 10*10*20 matrix
    pmf_zerodiag[[i]] <- sum_h1h2_zerodiag[[i]] / sum_h1h2rhok_zerodiag
  }
  result_list <- list("pmf_adj" = pmf_adj, 
                      "pmf_zerodiag" = pmf_zerodiag)
  return(result_list) 
}
