
source("~/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_yiting_gencor_CI/Code/journal_submission_GitHub_version/parametric_bootstrap_helper_functions.R")

######### function to get pmf #############
get_pmf <- function(h1_sq_hat, h2_sq_hat, rho_k_hat, adj_or_zerodiag, pmf_adj, pmf_zerodiag, bin_size){
  h1_sq_hat_index <- bin_index_of_h_sq_hat(h1_sq_hat, bin_size)
  h2_sq_hat_index <- bin_index_of_h_sq_hat(h2_sq_hat, bin_size)
  rho_k_hat_index <- bin_index_of_rho_hat(rho_k_hat, bin_size)
  my_pmf <- rep(0, 20)
  if(adj_or_zerodiag == 'adj'){
    #print('adj')
    pmf <- pmf_adj
  }else{
    #print('zerodiag')
    pmf <- pmf_zerodiag
  }
  for(i in 1:20){
    # pmf is i=1:20 given (h1_hat, h2_hat, rho_k_hat)
    my_pmf[i] = pmf[[i]][h1_sq_hat_index,h2_sq_hat_index,rho_k_hat_index]
  }
  return(my_pmf)
}

######### functions to get CI from empirical pmf #############
extend_right <- function(bins, rho_hat, percentile_hat, prob_from_hat, bin_size){
  cdf <- Reduce('+', bins, accumulate = TRUE)
  upper_CI_percentile <- percentile_hat + prob_from_hat
  upper_CI_index <- which(cdf >= upper_CI_percentile)[1]
  upper_CI_label <- index_to_label_rho(upper_CI_index, bin_size)
  high <- cdf[which(cdf >= upper_CI_percentile)[1]]
  low <- ifelse(which(cdf >= upper_CI_percentile)[1] == 1, 0, cdf[which(cdf >= upper_CI_percentile)[1]-1])
  upper_CI <- approx(x = c(low, high), y = c(upper_CI_label, upper_CI_label+bin_size), 
                     xout = upper_CI_percentile, method="linear")$y
  return(upper_CI)
}

extend_left <- function(bins, rho_hat, curr_bin_prob_left, prob_from_hat, bin_size){
  bin_hat_index <- bin_index_of_rho_hat(rho_hat, bin_size)
  if(bin_hat_index == 1){ # rho_hat is in the first bin
    low <- curr_bin_prob_left
    high <- 0
    lower_CI <- approx(x = c(low, high), y = c(-1, rho_hat), 
                       xout = prob_from_hat, method="linear")$y
  }else{
    bins_left <- bins[1:bin_hat_index-1]
    bins_left_cdf_rev <- c(c(0), Reduce('+', rev(bins_left), accumulate = TRUE)) # add a zero term for the left half of the current bin 
    cdf_rev_from_hat <- bins_left_cdf_rev + curr_bin_prob_left
    lower_CI_rev_index <- which(cdf_rev_from_hat >= prob_from_hat)[1]
    if(lower_CI_rev_index == 1){ # lower_CI is in the same bin as rho_hat
      low <- curr_bin_prob_left
      high <- 0
      bin_label <- index_to_label_rho(bin_hat_index, bin_size)
      lower_CI <- approx(x = c(low, high), y = c(bin_label, rho_hat), 
                         xout = prob_from_hat, method="linear")$y
    }else{
      low <- cdf_rev_from_hat[lower_CI_rev_index]
      high <- cdf_rev_from_hat[lower_CI_rev_index - 1]
      lower_CI_index <- bin_hat_index - lower_CI_rev_index + 1
      lower_CI_label <- index_to_label_rho(lower_CI_index, bin_size)
      lower_CI <- approx(x = c(low, high), y = c(lower_CI_label, lower_CI_label+bin_size), 
                         xout = prob_from_hat, method="linear")$y
    }
  }
  return(lower_CI)
}

get_CI <- function(bins, rho_hat, alpha, bin_size){ #bins = my_pmf, rho_hat = rho_k_hat, alpha = 0.05, bin_size = 0.1
  percentile_hat <- get_percentile_of_rho_hat(rho_hat, bin_size, bins)$percentile_hat
  curr_bin_prob_left <- get_percentile_of_rho_hat(rho_hat, bin_size, bins)$curr_bin_prob_left
  
  ### 3 cases of percentile: p < (1-alpha)/2, p > 1 - (1-alpha)/2, in between
  if(percentile_hat < (1-alpha)/2){
    ## 1. left edge: p < (1-alpha)/2
    # lower CI: the lower bound of the first bin with prob != 0
    # upper CI: extend (1-alpha-p) to the right
    prob_from_hat <- 1-alpha-percentile_hat
    lower_CI <- index_to_label_rho(which(bins>0)[1], bin_size) #as.numeric(names(which(bins>0)[1])) 
    upper_CI <- extend_right(bins, rho_hat, percentile_hat, prob_from_hat, bin_size)
    
  }else if(percentile_hat > 1 - (1-alpha)/2){
    ## 2. right edge: p > 1 - (1-alpha)/2
    # lower CI: extend (p - alpha) to the left
    # upper CI: the upper bound of the last bin with prob != 0
    prob_from_hat <- percentile_hat-alpha
    lower_CI <- extend_left(bins, rho_hat, curr_bin_prob_left, prob_from_hat, bin_size)
    upper_CI <- index_to_label_rho(rev(which(bins>0))[1], bin_size) + bin_size #as.numeric(names(rev(which(bins>0))[1])) + 0.01
    
  }else{
    ## 3. main case
    # lower CI: extend (1-alpha)/2 to the left
    # upper CI: extend (1-alpha)/2 to the right
    prob_from_hat <- (1-alpha)/2
    lower_CI <- extend_left(bins, rho_hat, curr_bin_prob_left, prob_from_hat, bin_size)
    upper_CI <- extend_right(bins, rho_hat, percentile_hat, prob_from_hat, bin_size)
  }
  CI_df <- data.frame(lower=lower_CI, upper=upper_CI)
  return(CI_df)
}

######### functions to get CI from beta #############
# helper functions
transform_to_beta_scale <- function(x){ (x + 1)/2}
transform_from_beta_scale <- function(x){ x*2 -1}
# get the gencor_pmf_values with actual count (big bins)
get_gencor_pmf_values <- function(my_pmf, n=10000){
  actual_counts <- round(n * my_pmf)
  sum(actual_counts) # not 10,000, but it's ok
  actual_values <- seq(-0.95, 0.95, by=0.1)
  gencor_pmf_values <- vector()
  for(i in 1:length(actual_values)){
    gencor_pmf_values <- append(gencor_pmf_values, rep(actual_values[i], actual_counts[i]))
  }
  return(gencor_pmf_values)
}

# return the CI of gencor [-1,1] based on beta distribution
get_CI_from_beta <- function(rho_est, shape1, shape2, coverage){
  # change rho_est into beta dist
  rho_est_beta <- transform_to_beta_scale(rho_est)
  # compute the probability in each site of val
  prob_left <- pbeta(rho_est_beta, shape1, shape2)
  prob_right <- 1-prob_left
  if (prob_left < coverage/2){ # case1: left edge
    left_point <- 0
    right_point <- qbeta(coverage, shape1, shape2)
  } else {
    if (prob_right < coverage/2){. # case2: right edge
      right_point <- 1
      left_point <- qbeta(1 - coverage, shape1, shape2)
    } else{ # case 3: we have sufficient "probabliity" on both sides
      left_point <- qbeta(prob_left - coverage/2 , shape1, shape2)
      right_point <- qbeta(coverage/2 + prob_left, shape1, shape2)
    }
  }
  CI_in_beta <- c(lower = left_point, upper = right_point)
  CI <- transform_from_beta_scale(CI_in_beta)
  return(CI)
}

######### function to get p-value from beta #############
# compute p-value (helper functions)
zero_in_CI <- function(CI){
  return(0 > CI[1] & 0 < CI[2])
}
zero_on_CI_tip <- function(CI){
  return(0 == CI[1] | 0 == CI[2])
}
binary_search_pval_beta <- function(rho_est, shape1, shape2, coverage_low, coverage_high, error_level){
  # the original coverage_low doesn't cover 0
  # the original coverage_high covers 0
  coverage_mid <- (coverage_low + coverage_high)/2
  if(coverage_high -coverage_low < error_level){
    return(1 - coverage_mid)
  }
  CI_mid <- get_CI_from_beta(rho_est, shape1, shape2, coverage_mid)
  if(zero_on_CI_tip(CI_mid)){
    return(1 - coverage_mid)
  }else if(zero_in_CI(CI_mid)){
    return(binary_search_pval_beta(rho_est, shape1, shape2, coverage_low, coverage_mid, error_level))
  }else{
    return(binary_search_pval_beta(rho_est, shape1, shape2, coverage_mid, coverage_high, error_level))
  }
}

# compute p-value
beta_pval <- function(rho_est, shape1, shape2, 
                      precision_cutoff=0.01, 
                      error_level_small_p=1e-12, 
                      error_level_large_p=0.01){
  # step 1: see p-value >=< 0.01
  CI <- get_CI_from_beta(rho_est, shape1, shape2, 1-precision_cutoff)
  if(zero_on_CI_tip(CI)){ # 0 is on the tip of 95% CI
    return(precision_cutoff)
  }else{
    if(zero_in_CI(CI)){ # 0 is covered in 95% CI => p-value > 0.01
      # do binary search
      coverage_low <- 0
      coverage_high <- 0.99
      return(binary_search_pval_beta(rho_est, shape1, shape2, coverage_low, coverage_high, error_level_large_p))
    }else{ # 0 is not covered in 95% CI => p-value < 0.01
      # do binary search
      coverage_low <- 0.99
      coverage_high <- 1
      return(binary_search_pval_beta(rho_est, shape1, shape2, coverage_low, coverage_high, error_level_small_p))
    }
  }
}

######### wrapper function to get CI from pmf, CI from beta, p-value from beta #############
wrapper_ftn <- function(h1_sq_hat, h2_sq_hat, rho_k_hat, adj_or_zerodiag, pmf_adj, pmf_zerodiag, bin_size=0.1, alpha=0.05, get_pval=T){
  # CI from pmf
  my_pmf <- get_pmf(h1_sq_hat,h2_sq_hat,rho_k_hat,adj_or_zerodiag, pmf_adj, pmf_zerodiag, bin_size)
  CI_from_pmf <- get_CI(my_pmf, rho_k_hat, alpha, bin_size) # get CI from pmf
  
  # CI from beta
  gencor_pmf_values <- get_gencor_pmf_values(my_pmf)
  
  if (length(unique(gencor_pmf_values)) == 1){
    cur_val <- unique(gencor_pmf_values)
    gencor_pmf_values <- c(rep(cur_val, length(gencor_pmf_values) - 20), 
                           rep(cur_val - bin_size/3, 10), 
                           rep(cur_val + bin_size/3, 10))
  }
  
  beta_dist <- ebeta(transform_to_beta_scale(gencor_pmf_values))
  shape1 <- beta_dist$parameters[["shape1"]]
  shape2 <- beta_dist$parameters[["shape2"]]
  CI_from_beta <- get_CI_from_beta(rho_k_hat, shape1, shape2, 1-alpha)

  # p-value from beta
  if(get_pval){
    pval_from_beta <- beta_pval(rho_k_hat, shape1, shape2)
  }else{
    pval_from_beta <- NA
  }

  result_list <- list("CI_from_pmf" = CI_from_pmf,
                      "CI_from_beta" = CI_from_beta,
                      "pval_from_beta" = pval_from_beta,
                      "pmf" = my_pmf,
                      "beta_shape1" = shape1,
                      "beta_shape2" = shape2)
  return(result_list)
  
  
   
}
