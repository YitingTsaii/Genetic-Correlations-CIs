library(psych)
library(mvtnorm)
library(hash)
library(boot)
library(Matrix) 
library(parallel) #to test how many CPUs are used
library(Rcpp)
library(mvtnorm)
library(mvnfast)
library(MASS)

library(peakRAM) #for measuring RAM usage 

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

#https://cran.r-project.org/web/packages/mvnfast/vignettes/mvnfast.html
# turn off BLAS parallelism for rmvn
library(RhpcBLASctl)
blas_set_num_threads(1) 

sourceCpp("gencor_functions_multithreading.cpp")
source("parametric_bootstrap_helper_functions.R")

args <- commandArgs(trailingOnly=TRUE)
job_idx <- as.numeric(args[1])
kinMat_block_path <- args[2] 
kin_mtrx_path <- args[3]
output_path <- args[4]
bin_size <- as.numeric(args[5])
rho_e <- as.numeric(args[6])
n.sim <- as.numeric(args[7]) 
num_blocks <- as.numeric(args[8]) 


#performance measuring files out
res_df_mvnfast_by_block_path <- args[9]
res_df_mvnfast_over_all_blocks_path <- args[10]

res_df_calcGC_KMZ_by_block_path <- args[11]
res_df_calcGC_KMA_by_block_path <- args[12]
res_df_calcHE_1_by_block_path <- args[13]
res_df_calcHE_2_by_block_path <- args[14]

res_df_calcGC_KMZ_final_over_all_blocks_path <- args[15]
res_df_calcGC_KMA_final_over_all_blocks_path <- args[16]
res_df_calcHE_1_final_over_all_blocks_path <- args[17]
res_df_calcHE_2_final_over_all_blocks_path <- args[18]



print("the number of cores detected inside R script is: ")
print(parallel::detectCores())

# measure overall time
print("started the clock: ")
ptm <- proc.time()

#load the kinship matrix block
RAM_file_load_1 <- peakRAM(kinMat_block <- readRDS(kinMat_block_path))
print(paste0("RAM usage RAM_file_load_1 is: ", RAM_file_load_1))

#load the kinship matrix
RAM_file_load_2 <- peakRAM(kinMat <- readRDS(kin_mtrx_path))
print(paste0("RAM usage RAM_file_load_2 is: ", RAM_file_load_2))

# bin size
start_rho <- -(1+(1-bin_size))/2
start_h <- bin_size/2
end <- (1+(1-bin_size))/2

# parameter settings
#rhos_k <- c(0) # for type 1 error analysis
rhos_k <- seq(start_rho, end, by=bin_size) # -0.95, -0.85, ..., -0.05, 0.05, ..., 0.85, 0.95 if bin_size = 0.1
h1s_sq <- seq(start_h, end, by=bin_size) # 0.05, 0.15, 0.25, …., 0.95
h2s_sq <- seq(start_h, end, by=bin_size) # 0.05, 0.15, 0.25, …., 0.95

ppl_n <- dim(kinMat)[1] 

paramspace<-as.vector(outer(rhos_k,h1s_sq,paste,sep=","))
paramspace<-as.vector(outer(paramspace,h2s_sq,paste,sep=","))

#length(paramspace) #2000 (when bin_size = 0.1)
#extract a triplet setting by job index
params <- paramspace[as.numeric(job_idx)] 
params <- as.numeric(unlist(strsplit(params, split=","))) 

# simulate data with known parameters (rho_k, h1_sq, h2_sq)
n <- ppl_n    
rho_k <- params[1]     
h1_sq <- params[2]     
h2_sq <- params[3]  

sprintf("number of simulations %d",n.sim)
sprintf("number of people in each simulation %d",n)
sprintf("rho_k %f",rho_k)
sprintf("h1 %f",h1_sq)
sprintf("h2 %f",h2_sq)

block_size <- dim(kinMat_block)[1] 
outcome_1 <- list()
outcome_2 <- list()

# construct covariance matrix between all outcomes (two outcomes for each individual),
# the end result (allmat.cov.all) is a 2nx2n matrix
allmat.cov.1 <- kinMat_block*h1_sq
allmat.cov.2 <- kinMat_block*h2_sq
allmat.cov.1.2 <- kinMat_block*sqrt(h1_sq)*sqrt(h2_sq)*rho_k
allmat.cov.all <- rbind(cbind(allmat.cov.1, allmat.cov.1.2),
                        cbind(allmat.cov.1.2, allmat.cov.2))


# construct covariance matrix for the independent errors 
# (i.e. not correlated between individuals)
indep.err.mat <- rbind(cbind((1-h1_sq)*diag(block_size), sqrt(1-h1_sq)*sqrt(1-h2_sq)*rho_e*diag(block_size)), 
                       cbind(sqrt(1-h1_sq)*sqrt(1-h2_sq)*rho_e*diag(block_size), (1-h2_sq)*diag(block_size))) 

#create empty df to store performance for mvnfast()
res_df_mvnfast <- data.frame(block_num = NA,
                             Function_Call = rep(NA, length(num_blocks)), 
                             Elapsed_Time_sec = NA, 
                             Total_RAM_Used_MiB = NA, 
                             Peak_RAM_Used_MiB = NA)

for (block in 1:num_blocks){
  
  combined_error_kinship_mat <- allmat.cov.all + indep.err.mat
  
  mu <- rep(0, 2*block_size)
  
  mcov <- combined_error_kinship_mat
  
  #measuring performance
  perf_measures <- peakRAM(b_plus_e_block <- mvnfast::rmvn(n.sim, mu, mcov, ncores = 8))

  #store performance results by block
  res_df_mvnfast[block,] <- data.frame(block_num = block,
                                       Function_Call = perf_measures$Function_Call,
                                       Elapsed_Time_sec = perf_measures$Elapsed_Time_sec,
                                       Total_RAM_Used_MiB = perf_measures$Total_RAM_Used_MiB,
                                       Peak_RAM_Used_MiB = perf_measures$Peak_RAM_Used_MiB)
  
  
  outcome_1[[block]] <- b_plus_e_block[1:n.sim, 1:block_size]
  
  outcome_2[[block]] <- b_plus_e_block[1:n.sim, (block_size + 1):(2*block_size)]
  
}

#save file
save(res_df_mvnfast, file = res_df_mvnfast_by_block_path)

#store_final_results (note res_df_mvnfast$Function_Call[1] means extract only first row for the name of the function so no repeats)
res_df_mvnfast_final <- data.frame(Function_Call = res_df_mvnfast$Function_Call[1], 
                                   Elapsed_Time_sec = Reduce('+', res_df_mvnfast$Elapsed_Time_sec), 
                                   Total_RAM_Used_MiB = Reduce('+', res_df_mvnfast$Total_RAM_Used_MiB), 
                                   Peak_RAM_Used_MiB = Reduce('+', res_df_mvnfast$Peak_RAM_Used_MiB))

#save file
save(res_df_mvnfast_final, file = res_df_mvnfast_over_all_blocks_path)

outcome_1_combined <- do.call(cbind, outcome_1)
outcome_2_combined <- do.call(cbind, outcome_2)

b_plus_e <- cbind(outcome_1_combined, outcome_2_combined)


# design matrix for the Haseman-Elston regression
Xs <- matrix(0,ncol = 2, nrow = n*n)
Xs[,2] <- c(kinMat) 
Xs[seq(1,n*n,n+1),1] <- 1
XsTXs <- t(Xs) %*% Xs #to optimize
rownames(XsTXs) = colnames(XsTXs) = c("error", "Kinship")

## XsTXs is a 2*2 matrix
XsTXsinv <- solve(XsTXs) 


# get kinMat.zero.diag and kinMat.adjusted0
# one does not use diagonal values, the other does.
# The word "adjusted" reflects that the Haseman-Elston regression
# accounts for estimation of the two variance components. 
kinMat.zero.diag <- kinMat; diag(kinMat.zero.diag) <- 0
wgts = XsTXsinv["Kinship",]
errorMat = diag(x = 1, nrow = nrow(kinMat), ncol = nrow(kinMat))
kinMat.adjusted0 = wgts["error"]*errorMat + wgts["Kinship"]*kinMat #to optimize


#trKK <- tr(kinMat.zero.diag %*% kinMat.zero.diag)#to optimize

#optimized version
kinMat_block_zero_diag <- kinMat_block
diag(kinMat_block_zero_diag) <- 0
trKK <- num_blocks*tr(kinMat_block_zero_diag %*% kinMat_block_zero_diag)



# vector to store estimates of genetic correlations
rho_k_hat_adj <- matrix(NA,nrow = n.sim, ncol = 1)
colnames(rho_k_hat_adj)<- "Kinship"
rho_k_hat_0 <- matrix(NA,nrow = n.sim, ncol = 1)
colnames(rho_k_hat_0)<- "Kinship"

# vectors to store estimates of h1_sq_hat and h2_sq_hat
h1_sq_hat <- rep(0, n.sim)
h2_sq_hat <- rep(0, n.sim)


#create empty df to store performance for calcGenCorMany_block_mtx_threads() using kinMat.zero.diag
res_df_calcGC_KMZ <- data.frame(block_num = NA,
                             Function_Call = rep(NA, length(num_blocks)), 
                             Elapsed_Time_sec = NA, 
                             Total_RAM_Used_MiB = NA, 
                             Peak_RAM_Used_MiB = NA)


#create empty df to store performance for calcGenCorMany_block_mtx_threads() using kinMat.adjusted0
res_df_calcGC_KMA <- data.frame(block_num = NA,
                              Function_Call = rep(NA, length(num_blocks)), 
                              Elapsed_Time_sec = NA, 
                              Total_RAM_Used_MiB = NA, 
                              Peak_RAM_Used_MiB = NA)


#create empty df to store performance for calcHeritabilityManyWThreads() using sVC.1.mat and kinMat.zero.diag
res_df_calcHE_1 <- data.frame(block_num = NA,
                              Function_Call = rep(NA, length(num_blocks)), 
                              Elapsed_Time_sec = NA, 
                              Total_RAM_Used_MiB = NA, 
                              Peak_RAM_Used_MiB = NA)


#create empty df to store performance for calcHeritabilityManyWThreads() using sVC.2.mat and kinMat.zero.diag
res_df_calcHE_2 <- data.frame(block_num = NA,
                              Function_Call = rep(NA, length(num_blocks)), 
                              Elapsed_Time_sec = NA, 
                              Total_RAM_Used_MiB = NA, 
                              Peak_RAM_Used_MiB = NA)

block_size <- 1000 

num_blocks <- n.sim/block_size #(10000/1000 = 10)

for (i in 1:num_blocks){
  
  
  block_start <- ((i-1)*block_size) + 1
  
  block_end <- i*block_size
  
  
  sVC.1.mat <- b_plus_e[block_start:block_end, 1:n]
  sVC.2.mat <- b_plus_e[block_start:block_end, (n+1):(2*n)]
  
  #measure performance and memory for each of the function calls
  perf_measures_1 <- peakRAM(gencor_zero_diag<- calcGenCorMany_block_mtx_threads(sVC.1.mat, sVC.2.mat, kinMat.zero.diag, n))
  

  #store performance results by block
  res_df_calcGC_KMZ[i,] <- data.frame(block_num = i,
                                       Function_Call = perf_measures_1$Function_Call,
                                       Elapsed_Time_sec = perf_measures_1$Elapsed_Time_sec,
                                       Total_RAM_Used_MiB = perf_measures_1$Total_RAM_Used_MiB,
                                       Peak_RAM_Used_MiB = perf_measures_1$Peak_RAM_Used_MiB)
  
  
  perf_measures_2 <- peakRAM(gencor_kinship_adjusted <- calcGenCorMany_block_mtx_threads(sVC.1.mat, sVC.2.mat, kinMat.adjusted0, n))

  #store performance results by block
  res_df_calcGC_KMA[i,] <- data.frame(block_num = i,
                                      Function_Call = perf_measures_2$Function_Call,
                                      Elapsed_Time_sec = perf_measures_2$Elapsed_Time_sec,
                                      Total_RAM_Used_MiB = perf_measures_2$Total_RAM_Used_MiB,
                                      Peak_RAM_Used_MiB = perf_measures_2$Peak_RAM_Used_MiB)
  
  
  perf_measures_3 <- peakRAM(heritabilities_h1_sq_hat <- calcHeritabilityManyWThreads(sVC.1.mat, kinMat.zero.diag))
  
  #store performance results by block
  res_df_calcHE_1[i,] <- data.frame(block_num = i,
                                      Function_Call = perf_measures_3$Function_Call,
                                      Elapsed_Time_sec = perf_measures_3$Elapsed_Time_sec,
                                      Total_RAM_Used_MiB = perf_measures_3$Total_RAM_Used_MiB,
                                      Peak_RAM_Used_MiB = perf_measures_3$Peak_RAM_Used_MiB)
  
  
  perf_measures_4 <- peakRAM(heritabilities_h2_sq_hat <- calcHeritabilityManyWThreads(sVC.2.mat, kinMat.zero.diag))

  #store performance results by block
  res_df_calcHE_2[i,] <- data.frame(block_num = i,
                                    Function_Call = perf_measures_4$Function_Call,
                                    Elapsed_Time_sec = perf_measures_4$Elapsed_Time_sec,
                                    Total_RAM_Used_MiB = perf_measures_4$Total_RAM_Used_MiB,
                                    Peak_RAM_Used_MiB = perf_measures_4$Peak_RAM_Used_MiB)
  
  
  
  rho_k_hat_0[block_start:block_end, ] <- gencor_zero_diag
  rho_k_hat_adj[block_start:block_end,] <- gencor_kinship_adjusted
  h1_sq_hat[block_start:block_end] <- heritabilities_h1_sq_hat
  h2_sq_hat[block_start:block_end] <- heritabilities_h2_sq_hat
  
}

#save results by blocks 
save(res_df_calcGC_KMZ, file = res_df_calcGC_KMZ_by_block_path)
save(res_df_calcGC_KMA, file = res_df_calcGC_KMA_by_block_path)
save(res_df_calcHE_1, file = res_df_calcHE_1_by_block_path)
save(res_df_calcHE_2, file = res_df_calcHE_2_by_block_path)


#store_final_results (note Function_Call[1] means extract only first raw for the name of the function so no repeats)
res_df_calcGC_KMZ_final <- data.frame(Function_Call = res_df_calcGC_KMZ$Function_Call[1], 
                                      Elapsed_Time_sec = Reduce('+', res_df_calcGC_KMZ$Elapsed_Time_sec), 
                                      Total_RAM_Used_MiB = Reduce('+', res_df_calcGC_KMZ$Total_RAM_Used_MiB), 
                                      Peak_RAM_Used_MiB = Reduce('+', res_df_calcGC_KMZ$Peak_RAM_Used_MiB))

res_df_calcGC_KMA_final <- data.frame(Function_Call = res_df_calcGC_KMA$Function_Call[1], 
                                      Elapsed_Time_sec = Reduce('+', res_df_calcGC_KMA$Elapsed_Time_sec), 
                                      Total_RAM_Used_MiB = Reduce('+', res_df_calcGC_KMA$Total_RAM_Used_MiB), 
                                      Peak_RAM_Used_MiB = Reduce('+', res_df_calcGC_KMA$Peak_RAM_Used_MiB))

res_df_calcHE_1_final <- data.frame(Function_Call = res_df_calcHE_1$Function_Call[1], 
                                      Elapsed_Time_sec = Reduce('+', res_df_calcHE_1$Elapsed_Time_sec), 
                                      Total_RAM_Used_MiB = Reduce('+', res_df_calcHE_1$Total_RAM_Used_MiB), 
                                      Peak_RAM_Used_MiB = Reduce('+', res_df_calcHE_1$Peak_RAM_Used_MiB))

res_df_calcHE_2_final <- data.frame(Function_Call = res_df_calcHE_2$Function_Call[1], 
                                      Elapsed_Time_sec = Reduce('+', res_df_calcHE_2$Elapsed_Time_sec), 
                                      Total_RAM_Used_MiB = Reduce('+', res_df_calcHE_2$Total_RAM_Used_MiB), 
                                      Peak_RAM_Used_MiB = Reduce('+', res_df_calcHE_2$Peak_RAM_Used_MiB))

#save final results (over all blocks) 
save(res_df_calcGC_KMZ_final, file = res_df_calcGC_KMZ_final_over_all_blocks_path)
save(res_df_calcGC_KMA_final, file = res_df_calcGC_KMA_final_over_all_blocks_path)
save(res_df_calcHE_1_final, file = res_df_calcHE_1_final_over_all_blocks_path)
save(res_df_calcHE_2_final, file = res_df_calcHE_2_final_over_all_blocks_path)



# output file function
save_file <- function(rho_k_hat_adj, rho_k_hat_0, h1_sq_hat, h2_sq_hat, rho_k, h1_sq, h2_sq){
  outfname=sprintf(paste0(output_path, 'sim_genecorr_run-%d_nsim-%d_rhok-%.3f_h1sq-%.3f_h2sq-%.3f.RData'),
                   job_idx,n.sim,rho_k,h1_sq,h2_sq)
  
  # calculate 3D bins
  bins_3D_adj <- cal_3D_bins(n.sim, h1_sq_hat, h2_sq_hat, rho_k_hat_adj, bin_size)
  bins_3D_zerodiag <- cal_3D_bins(n.sim, h1_sq_hat, h2_sq_hat, rho_k_hat_0, bin_size)
  
  # create rres list
  rres <- list()
  rres$rho_k_hat_adj <- rho_k_hat_adj
  rres$rho_k_hat_0 <- rho_k_hat_0
  rres$h1_sq_hat <- h1_sq_hat
  rres$h2_sq_hat <- h2_sq_hat
  rres$bins_3D_adj <- bins_3D_adj
  rres$bins_3D_zerodiag <- bins_3D_zerodiag
  rres$rho_k <- rho_k
  rres$h1_sq <- h1_sq
  rres$h2_sq <- h2_sq
  
  save(rres,file=outfname)
}


save_file(rho_k_hat_adj, rho_k_hat_0, h1_sq_hat, h2_sq_hat, rho_k, h1_sq, h2_sq)
print("Done.")


end_time <- proc.time() - ptm
print(paste0("The overall runtime is: ", end_time))

