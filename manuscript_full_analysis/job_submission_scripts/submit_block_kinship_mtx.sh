#!/bin/bash

kin_mtrx_path="/Data/Kinship_matrix/freeze.8.JHS_mtx/JHS_kinship_matrix_from_dense_TOPMed_freeze_8.rds"
output_path="/Data/Kinship_matrix/freeze.8.JHS_mtx_3X/JHS_kinship_matrix_from_dense_TOPMed_freeze_8_3X.rds" 
mult_coef=3 # 3X or 3X


Rscript 2_block_kinship.R $kin_mtrx_path $output_path $mult_coef