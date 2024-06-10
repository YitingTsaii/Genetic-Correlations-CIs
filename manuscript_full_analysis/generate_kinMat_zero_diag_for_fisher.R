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

args <- commandArgs(trailingOnly=TRUE)

kin_mtrx_path <- args[1]
output_path <- args[2]

#load the kinship matrix
kinMat <- readRDS(kin_mtrx_path)

kinMat.zero.diag <- kinMat; 

diag(kinMat.zero.diag) <- 0


saveRDS(kinMat.zero.diag, output_path)
