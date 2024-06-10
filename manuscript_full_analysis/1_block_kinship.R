library(Matrix) 



args <- commandArgs(trailingOnly=TRUE)
kin_mtrx_path <- args[1]
output_path <- args[2] 
mult_coef <- args[3] 



kinMat <- readRDS(kin_mtrx_path)

new_kin_mtx <- kronecker(diag(mult_coef), kinMat)


colnames(new_kin_mtx) <- paste0("P", c(1:ncol(new_kin_mtx)))

rownames(new_kin_mtx) <- paste0("P", c(1:ncol(new_kin_mtx)))

saveRDS(new_kin_mtx, output_path)




