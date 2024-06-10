library(psych)
library(mvtnorm)
library(hash)
library(boot)
library(Matrix) 
library(parallel) 
library(Rcpp)
library(mvtnorm)
library(mvnfast)
library(MASS)

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

#https://cran.r-project.org/web/packages/mvnfast/vignettes/mvnfast.html
# We also need to turn off BLAS parallelism for running rmvn
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
n.sim <- as.numeric(args[7]) # number of simulations = N
num_blocks <- as.numeric(args[8]) 


print("the number of cores detected inside R script is: ")
print(parallel::detectCores())

#load the kinship matrix block
kinMat_block <- readRDS(kinMat_block_path)


#load the kinship matrix
kinMat <- readRDS(kin_mtrx_path)


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


for (block in 1:num_blocks){
  
  combined_error_kinship_mat <- allmat.cov.all + indep.err.mat
  
  mu <- rep(0, 2*block_size)
  
  mcov <- combined_error_kinship_mat
  
  b_plus_e_block <- mvnfast::rmvn(n.sim, mu, mcov, ncores = 8)
  
  outcome_1[[block]] <- b_plus_e_block[1:n.sim, 1:block_size]
  
  outcome_2[[block]] <- b_plus_e_block[1:n.sim, (block_size + 1):(2*block_size)]
  
}


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
XsTXsinv <- solve(XsTXs) #to optimize (most likely ok, double check)


# get kinMat.zero.diag and kinMat.adjusted0
# one does not use diagonal values, the other does.
# The word "adjusted" reflects that the Haseman-Elston regression
# accounts for estimation of the two variance components. 
kinMat.zero.diag <- kinMat; diag(kinMat.zero.diag) <- 0
wgts = XsTXsinv["Kinship",]
errorMat = diag(x = 1, nrow = nrow(kinMat), ncol = nrow(kinMat))
kinMat.adjusted0 = wgts["error"]*errorMat + wgts["Kinship"]*kinMat #to optimize


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


#we hardwrite this value just because we know that we have a total of 
#10000 simulations and we want to split them into 10 round of a 1000
block_size <- 1000

num_blocks <- n.sim/block_size #(10000/1000 = 10)



for (i in 1:num_blocks){
  
  
  block_start <- ((i-1)*block_size) + 1
  
  block_end <- i*block_size
  
  
  sVC.1.mat <- b_plus_e[block_start:block_end, 1:n]
  sVC.2.mat <- b_plus_e[block_start:block_end, (n+1):(2*n)]
  
  
  gencor_zero_diag<- calcGenCorMany_block_mtx_threads(sVC.1.mat, sVC.2.mat, kinMat.zero.diag, n)
  gencor_kinship_adjusted <- calcGenCorMany_block_mtx_threads(sVC.1.mat, sVC.2.mat, kinMat.adjusted0, n)
  
  heritabilities_h1_sq_hat <- calcHeritabilityManyWThreads(sVC.1.mat, kinMat.zero.diag)
  heritabilities_h2_sq_hat <- calcHeritabilityManyWThreads(sVC.2.mat, kinMat.zero.diag)
  
  rho_k_hat_0[block_start:block_end, ] <- gencor_zero_diag
  rho_k_hat_adj[block_start:block_end,] <- gencor_kinship_adjusted
  h1_sq_hat[block_start:block_end] <- heritabilities_h1_sq_hat
  h2_sq_hat[block_start:block_end] <- heritabilities_h2_sq_hat
  
}



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