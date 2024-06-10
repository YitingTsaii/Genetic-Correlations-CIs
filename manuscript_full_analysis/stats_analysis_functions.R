
source("~/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_yiting_gencor_CI/Code/journal_submission_GitHub_version/parametric_bootstrap_helper_functions.R")
source("~/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/2022_yiting_gencor_CI/Code/journal_submission_GitHub_version/compute_CI_functions.R")

######### function to get cover_prob_all, CI_width_all for pmf #############
get_summary_df_pmf <- function(path, pmf_adj, pmf_zerodiag){
  # path to the folder with 2000 simulation results. For example.
  # path <- "/data/linkage/JHS/Projects/2022_bootstrap_gencor_CIs/Simulations/0221104_new_kinship_C" 
  
  #sanity check
  print(pmf_adj)
  print(pmf_zerodiag) 
  
  # average over all the 100 (h1,h2) settings
  cover_prob_all <- data.frame(matrix(seq(-0.95, 0.95, by=0.1)))
  CI_width_all <- data.frame(matrix(seq(-0.95, 0.95, by=0.1)))
  colnames(cover_prob_all) <- "true_rho_k"
  colnames(CI_width_all) <- "true_rho_k"
  
  #cover_prob_vec <- vector()
  for(true_h1_sq in seq(0.05,0.95,by=0.1)){
    for(true_h2_sq in seq(0.05,0.95,by=0.1)){
      print(paste0("true_h1_sq=",true_h1_sq))
      print(paste0("true_h2_sq=",true_h2_sq))
      cover_prob_mat <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(cover_prob_mat) <- c("true_rho_k", "cover_prob", "ave_width_CI", "ave_lower_CI", "ave_upper_CI")
      i <- 1
      
      for(true_rho_k in seq(-0.95, 0.95, 0.1)){
        print(paste0("true_rho_k=",true_rho_k))
        pattern <- paste0("rhok-", round(true_rho_k,3),"0_h1sq-",round(true_h1_sq,2),"0_h2sq-",round(true_h2_sq,2),"0") 
        
        file_name <- list.files(path=path, pattern = pattern)

    
        load(paste0(path, "/", file_name))
        
        CI_matrix <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(CI_matrix) <- c("rho_hat","lower_CI", "upper_CI")
        for(j in 1:10000){
          # 10,000 rho_hat in total for each true_rho_k
          h1_sq_hat <- rres$h1_sq_hat[j]
          h2_sq_hat <- rres$h2_sq_hat[j]
          rho_k_hat <- rres$rho_k_hat_adj[j]
          
          if(is.na(rho_k_hat)){  # skip the iteration if rho_k_hat is NA
            next
          }else if(rho_k_hat < -1){
            rho_k_hat <- -1
          }else if(rho_k_hat > 1){
            rho_k_hat <- 1
          }
          
          wrapper_result <- wrapper_ftn(h1_sq_hat, h2_sq_hat, rho_k_hat, 'adj', pmf_adj, pmf_zerodiag, bin_size=0.1, alpha=0.05, get_pval=F)
          CI <- wrapper_result$CI_from_pmf
          
          new_row <- t(matrix(c(rho_k_hat,CI$lower,CI$upper)))
          colnames(new_row) <- c("rho_hat","lower_CI", "upper_CI")
          CI_matrix <- rbind(CI_matrix, new_row)
        }
        
        CI_matrix$CI_width = CI_matrix$upper_CI - CI_matrix$lower_CI
        CI_summary <- apply(CI_matrix, 2, mean)
        in_CI_result <- apply(CI_matrix, 1, function(x) true_rho_in_CI(true_rho_k, x[2], x[3]))
        cover_prob <- sum(in_CI_result)/length(in_CI_result)
        #cover_prob_vec <- append(cover_prob_vec, cover_prob)
        newdata <- t(matrix(c(true_rho_k,cover_prob,CI_summary[4],CI_summary[2],CI_summary[3])))
        colnames(newdata) <- c("true_rho_k", "cover_prob", "ave_width_CI", "ave_lower_CI", "ave_upper_CI")
        cover_prob_mat <- rbind(cover_prob_mat, newdata)
        #
        print(i)
        i <- i+1
      }
      
      # cbind the coverage prob
      new_col <- data.frame(matrix(cover_prob_mat$cover_prob))
      colnames(new_col) <- paste0("h1_",true_h1_sq,"_h2_",true_h2_sq)
      cover_prob_all <- cbind(cover_prob_all, new_col)
      
      # cbind the average CI width
      new_col_width <- data.frame(matrix(cover_prob_mat$ave_width_CI))
      colnames(new_col_width) <- paste0("h1_",true_h1_sq,"_h2_",true_h2_sq)
      CI_width_all <- cbind(CI_width_all, new_col_width)
    } 
  }
  result_list <- list("cover_prob_all" = cover_prob_all, 
                      "CI_width_all" = CI_width_all)
  return(result_list) 
}



######### function to get cover_prob_all, CI_width_all for beta #############
get_summary_df_beta <- function(path, pmf_adj, pmf_zerodiag){
  # path to the folder with 2000 simulation results. For example.
  # path <- "/data/linkage/JHS/Projects/2022_bootstrap_gencor_CIs/Simulations/0221104_new_kinship_C" 
  cover_prob_all <- data.frame(matrix(seq(-0.95, 0.95, by=0.1)))
  CI_width_all <- data.frame(matrix(seq(-0.95, 0.95, by=0.1)))
  colnames(cover_prob_all) <- "true_rho_k"
  colnames(CI_width_all) <- "true_rho_k"
  
  for(true_h1_sq in seq(0.05,0.95,by=0.1)){
    for(true_h2_sq in seq(0.05,0.95,by=0.1)){
      #print(paste0("true_h1_sq=",true_h1_sq))
      #print(paste0("true_h2_sq=",true_h2_sq))
      cover_prob_mat <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(cover_prob_mat) <- c("true_rho_k", "cover_prob", "ave_width_CI", "ave_lower_CI", "ave_upper_CI")
      i <- 1
      
      for(true_rho_k in seq(-0.95, 0.95, 0.1)){
        print(paste0("true_rho_k=",true_rho_k))
        pattern <- paste0("rhok-", round(true_rho_k,3),"0_h1sq-",round(true_h1_sq,2),"0_h2sq-",round(true_h2_sq,2),"0") 
        file_name <- list.files(path=path, pattern = pattern)
        load(paste0(path, "/", file_name))

        
        CI_matrix <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(CI_matrix) <- c("rho_hat","lower_CI", "upper_CI")
        for(j in 1:10000){
          # 10,000 rho_hat in total for each true_rho_k
          h1_sq_hat <- rres$h1_sq_hat[j]
          h2_sq_hat <- rres$h2_sq_hat[j]
          rho_k_hat <- rres$rho_k_hat_adj[j]
          
          if(is.na(rho_k_hat)){  # skip the iteration if rho_k_hat is NA
            next
          }else if(rho_k_hat < -1){
            rho_k_hat <- -1
          }else if(rho_k_hat > 1){
            rho_k_hat <- 1
          }
          
          # get_CI
          alpha <- 0.05
          wrapper_result <- wrapper_ftn(h1_sq_hat, h2_sq_hat, rho_k_hat, 'adj', pmf_adj, pmf_zerodiag, bin_size=0.1, alpha=0.05, get_pval=F)
          CI_beta <- wrapper_result$CI_from_beta
          
          new_row <- t(matrix(c(rho_k_hat,CI_beta['lower'],CI_beta['upper'])))
          colnames(new_row) <- c("rho_hat","lower_CI", "upper_CI")
          CI_matrix <- rbind(CI_matrix, new_row)
        }
        
        CI_matrix$CI_width = CI_matrix$upper_CI - CI_matrix$lower_CI
        CI_summary <- apply(CI_matrix, 2, mean)
        
        
        in_CI_result <- apply(CI_matrix, 1, function(x) true_rho_in_CI(true_rho_k, x[2], x[3]))
        cover_prob <- sum(in_CI_result)/length(in_CI_result)
        newdata <- t(matrix(c(true_rho_k,cover_prob,CI_summary[4],CI_summary[2],CI_summary[3])))
        colnames(newdata) <- c("true_rho_k", "cover_prob", "ave_width_CI", "ave_lower_CI", "ave_upper_CI")
        cover_prob_mat <- rbind(cover_prob_mat, newdata)
        #
        print(i)
        i <- i+1
      }
      
      # cbind the coverage prob
      new_col <- data.frame(matrix(cover_prob_mat$cover_prob))
      colnames(new_col) <- paste0("h1_",true_h1_sq,"_h2_",true_h2_sq)
      cover_prob_all <- cbind(cover_prob_all, new_col)
      
      # cbind the average CI width
      new_col_width <- data.frame(matrix(cover_prob_mat$ave_width_CI))
      colnames(new_col_width) <- paste0("h1_",true_h1_sq,"_h2_",true_h2_sq)
      CI_width_all <- cbind(CI_width_all, new_col_width)
    } 
  }
  result_list <- list("cover_prob_all" = cover_prob_all, 
                      "CI_width_all" = CI_width_all)
  return(result_list) 
}


######### function to get cover_prob_all, CI_width_all for fisher #############
fisher_transform <- function(x){
  if (x <= -1) x <- -0.9999
  if (x >= 1) x <- 0.9999
  0.5*log((1+x)/(1-x))
}

fisher_inv_transform <- function(x){
  (exp(2*x) -1)/(exp(2*x) + 1)
}

calc.fisher <- function(gencor, n, conf = 0.05){
  if(is.na(gencor)){
    return(list(low = NA, high = NA))
  }
  if((gencor < -1)|(gencor > 1)){
    return(list(low = NA, high = NA))
  }
  qnorm_conf <- qnorm(conf/2, lower.tail = FALSE )
  z <- fisher_transform(gencor)
  low.z <- z - qnorm_conf/sqrt(n-3)
  high.z <- z + qnorm_conf/sqrt(n-3)
  
  low = fisher_inv_transform(low.z); high = fisher_inv_transform(high.z)
  return(list(low=low,high=high))
}

get_summary_df_fisher <- function(path, kinMat.adj_n){
  # path to the folder with 2000 simulation results. For example.
  # path <- "/data/linkage/JHS/Projects/2022_bootstrap_gencor_CIs/Simulations/0221104_new_kinship_C" 
  
  # average over all the 100 (h1,h2) settings
  cover_prob_all <- data.frame(matrix(seq(-0.95, 0.95, by=0.1)))
  CI_width_all <- data.frame(matrix(seq(-0.95, 0.95, by=0.1)))
  colnames(cover_prob_all) <- "true_rho_k"
  colnames(CI_width_all) <- "true_rho_k"
  
  #cover_prob_vec <- vector()
  for(true_h1_sq in seq(0.05,0.95,by=0.1)){
    for(true_h2_sq in seq(0.05,0.95,by=0.1)){
      #print(paste0("true_h1_sq=",true_h1_sq))
      #print(paste0("true_h2_sq=",true_h2_sq))
      cover_prob_mat <- data.frame(matrix(ncol = 5, nrow = 0))
      colnames(cover_prob_mat) <- c("true_rho_k", "cover_prob", "ave_width_CI", "ave_lower_CI", "ave_upper_CI")
      i <- 1
      
      for(true_rho_k in seq(-0.95, 0.95, 0.1)){
        print(paste0("true_rho_k=",true_rho_k))
        pattern <- paste0("rhok-", round(true_rho_k,3),"0_h1sq-",round(true_h1_sq,2),"0_h2sq-",round(true_h2_sq,2),"0") 
        
        
        file_name <- list.files(path=path, pattern = pattern)

        load(paste0(path, "/", file_name))
        
        CI_matrix <- data.frame(matrix(ncol = 3, nrow = 0))
        colnames(CI_matrix) <- c("rho_hat","lower_CI", "upper_CI")
        #for(rho_hat in rres$rho_k_hat_adj){
        for(j in 1:10000){
          rho_k_hat <- rres$rho_k_hat_adj[j]
          if(is.na(rho_k_hat)){  # skip the iteration if rho_k_hat is NA
            next
          }
          CI <- unlist(calc.fisher(rho_k_hat, kinMat.adj_n, conf = 0.05))
          new_row <- t(matrix(c(rho_k_hat,CI[1],CI[2])))
          colnames(new_row) <- c("rho_hat","lower_CI", "upper_CI")
          CI_matrix <- rbind(CI_matrix, new_row) # will have some NAs due to the function of getting fisher CI (if gencor > 1 or < -1, return CI = [NA, NA])
        }
        # for fisher CI
        # take out the rows with NAs
        # 
        CI_matrix <- na.omit(CI_matrix)
        
        CI_matrix$CI_width = CI_matrix$upper_CI - CI_matrix$lower_CI
        CI_summary <- apply(CI_matrix, 2, mean)
        in_CI_result <- apply(CI_matrix, 1, function(x) true_rho_in_CI(true_rho_k, x[2], x[3]))
        cover_prob <- sum(in_CI_result)/length(in_CI_result)
        newdata <- t(matrix(c(true_rho_k,cover_prob,CI_summary[4],CI_summary[2],CI_summary[3])))
        colnames(newdata) <- c("true_rho_k", "cover_prob", "ave_width_CI", "ave_lower_CI", "ave_upper_CI")
        cover_prob_mat <- rbind(cover_prob_mat, newdata)
        #
        print(i)
        i <- i+1
      }
      
      # cbind the coverage prob
      new_col <- data.frame(matrix(cover_prob_mat$cover_prob))
      colnames(new_col) <- paste0("h1_",true_h1_sq,"_h2_",true_h2_sq)
      cover_prob_all <- cbind(cover_prob_all, new_col)
      
      # cbind the average CI width
      new_col_width <- data.frame(matrix(cover_prob_mat$ave_width_CI))
      colnames(new_col_width) <- paste0("h1_",true_h1_sq,"_h2_",true_h2_sq)
      CI_width_all <- cbind(CI_width_all, new_col_width)
    } 
  }
  result_list <- list("cover_prob_all" = cover_prob_all, 
                      "CI_width_all" = CI_width_all)
  return(result_list) 
}


true_rho_in_CI <- function(true_rho_k, lower_CI, upper_CI){
  if(true_rho_k >= lower_CI & true_rho_k <= upper_CI){
    return(TRUE)
  }else{
    return(FALSE)
  }
}