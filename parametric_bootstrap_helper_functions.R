######### functions to get rho & h labels and indices #############
# rho bins:  [-1,-0.9),[-0.9,-0.8),...,[0.8,0.9),[0.9,1]
# rho label:        -1,       -0.9,...,      0.8,    0.9
# rho index:         1,          2,...,       19,     20
label_to_index_rho <- function(bin_label, bin_size){
  mult_fact <- 1/bin_size
  bin_index <- mult_fact*(bin_label + 1) + 1
  return(round(bin_index)) # round bin_index to the nearest integer to avoid the floating number input of bin_label
} 

index_to_label_rho <- function(bin_index, bin_size){
  mult_fact <- 1/bin_size
  bin_label <- -1 + (bin_index-1) * bin_size 
  bin_label_round <- round(bin_label*mult_fact)/mult_fact # round to the nearest "bin_size" to avoid floating number
  return(bin_label_round)
}

bin_index_of_rho_hat <- function(rho_hat, bin_size){ 
  # inputs can be any real value, so we now add the illegal inputs to the edge bins
  if(rho_hat >= 1){  # including all the illegal estimates >= 1
    bin_hat_label <- 1-bin_size
  }else if(rho_hat < -1){  # including all the illegal estimates < -1
    bin_hat_label <- -1
  }else{
    mult_fact <- 1/bin_size
    bin_hat_label <- floor(rho_hat*mult_fact)/mult_fact # take the floor to the nearest 0.01 or 0.1 or 0.05 ... ("biz_size")
  }
  bin_hat_index <- label_to_index_rho(bin_hat_label, bin_size)
  return(bin_hat_index)
}

# for h1 and h2
# h bins:  [0,0.1),[0.1,0.2),...,[0.8,0.9),[0.9,1]
# h label:       0,      0.1,...,      0.8,    0.9
# h index:       1,        2,...,        9,     10
label_to_index_h <- function(bin_label, bin_size){
  mult_fact <- 1/bin_size
  bin_index <- mult_fact*bin_label + 1
  return(round(bin_index))
} 

index_to_label_h <- function(bin_index, bin_size){
  mult_fact <- 1/bin_size
  bin_label <- (bin_index - 1)*bin_size
  bin_label_round <- round(bin_label*mult_fact)/mult_fact # round to the nearest "bin_size" to avoid floating number
  return(bin_label_round)
}

bin_index_of_h_sq_hat <- function(h_sq_hat, bin_size){ 
  # inputs can be any real value, add the illegal inputs to the edge bins
  if(h_sq_hat >= 1){  # including all the illegal estimates >= 1
    bin_hat_label <- 1-bin_size
  }else if(h_sq_hat < 0){   #including all the illegal estimates < 0
    bin_hat_label <- 0
  }else{
    mult_fact <- 1/bin_size
    bin_hat_label <- floor(h_sq_hat*mult_fact)/mult_fact # round to the nearest 0.01 or 0.1 or 0.05 ... ("biz_size")
  }
  bin_hat_index <- label_to_index_h(bin_hat_label, bin_size)
  return(bin_hat_index)
}

# calculate 3D bins [h1_sq, h2_sq, rho_k]
cal_3D_bins <- function(n.sim, h1_sq_hat, h2_sq_hat, rho_k_hat, bin_size){
  bins_3D <- array(rep(0,10*10*20), dim=c(10,10,20))  #bins_3D[h1_sq, h2_sq, rho_k]
  for(i in 1:n.sim){
    if(is.na(rho_k_hat[i])){  # skip the iteration if rho_k_hat is NA
      next
    }else{
      x <- bin_index_of_h_sq_hat(h1_sq_hat[i], bin_size)
      y <- bin_index_of_h_sq_hat(h2_sq_hat[i], bin_size)
      z <- bin_index_of_rho_hat(rho_k_hat[i], bin_size)
      bins_3D[x,y,z] <- bins_3D[x,y,z] + 1
    }
  }
  return(bins_3D)
}

get_percentile_of_rho_hat <- function(rho_hat, bin_size, bins){
  ## 1. identify the bin 
  bin_hat_index <- bin_index_of_rho_hat(rho_hat, bin_size)
  bin_hat_label <- index_to_label_rho(bin_hat_index, bin_size)
  
  ## 2. interpolation
  acc_prob_lower <- ifelse(bin_hat_index == 1, 0, sum(bins[1:(bin_hat_index-1)]))
  acc_prob_upper <- sum(bins[1:(bin_hat_index)])
  
  if (rho_hat < -1){
    percentile_hat <- approx(x = c(bin_hat_label, bin_hat_label+bin_size), 
                             y = c(acc_prob_lower, acc_prob_upper), 
                             xout = -1, method="linear")$y
  } else if (rho_hat > 1){
    percentile_hat <- approx(x = c(bin_hat_label, bin_hat_label+bin_size), 
                             y = c(acc_prob_lower, acc_prob_upper), 
                             xout = 1, method="linear")$y
  } else {
    percentile_hat <- approx(x = c(bin_hat_label, bin_hat_label+bin_size), 
                             y = c(acc_prob_lower, acc_prob_upper), 
                             xout = rho_hat, method="linear")$y
  }
  
  curr_bin_prob_left <- percentile_hat - acc_prob_lower
  
  my_list <- list("percentile_hat" = percentile_hat, "curr_bin_prob_left" = curr_bin_prob_left)
  return(my_list) 
}


