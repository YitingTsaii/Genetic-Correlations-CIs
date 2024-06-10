#include <Rcpp.h>
#include <omp.h>
#include <mutex>
#include <vector>
#include <iostream>
#include <cmath>
#include <thread>
#include <utility> 
#include <functional> // Include for std::function
#include <iomanip> // For std::fixed and std::setprecision
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::export]]
// Function to calculate genetic correlations based on multiple simulations at once using threads



Eigen::VectorXd calcGenCorMany_block_mtx_threads(const Eigen::MatrixXd& sVC1Mat,
                                                 const Eigen::MatrixXd& sVC2Mat,
                                                 const Eigen::MatrixXd& kinMat,
                                                 const int blockSize) {
  

  
  const int numBlocks = kinMat.cols() / blockSize; 
  
  
  const unsigned int numThreads = std::min(numBlocks, static_cast<int>(std::thread::hardware_concurrency()));
  
  
  std::cout << "the number of threads available in C++ is: " << std::thread::hardware_concurrency() << std::endl;
  
  std::vector<Eigen::VectorXd> arg1(numBlocks); 
  std::vector<Eigen::VectorXd> arg2(numBlocks);
  std::vector<Eigen::VectorXd> arg3(numBlocks);
  
  
  std::function<void(int)> computeBlock = [&](int blockIndex) {
    
    
    
    Eigen::MatrixXd kinMat_block = kinMat.block(blockIndex * blockSize, blockIndex * blockSize, blockSize, blockSize); 
    
    //only extract a block (of blockSize)
    Eigen::MatrixXd sVC1Mat_block = sVC1Mat.block(0, blockIndex * blockSize, sVC1Mat.rows(), blockSize);
    
    Eigen::MatrixXd sVC2Mat_block = sVC2Mat.block(0, blockIndex * blockSize, sVC2Mat.rows(), blockSize);
    
    Eigen::MatrixXd kaVC1MatBlock = sVC1Mat_block * kinMat_block;
    
    Eigen::MatrixXd kaVC2MatBlock = sVC2Mat_block * kinMat_block; 
    
    
    Eigen::MatrixXd chnk1 = kaVC1MatBlock * sVC2Mat_block.transpose();
    Eigen::MatrixXd chnk2 = kaVC1MatBlock * sVC1Mat_block.transpose();
    Eigen::MatrixXd chnk3 = kaVC2MatBlock * sVC2Mat_block.transpose();
    
    arg1[blockIndex] = chnk1.diagonal();
    arg2[blockIndex] = chnk2.diagonal();
    arg3[blockIndex] = chnk3.diagonal();
    
    
  };
  
  
  std::vector<std::thread> threads(numThreads);
  
  for (int i = 0; i < numBlocks; i += numThreads) {
    for (unsigned int j = 0; j < numThreads && i + j < numBlocks; ++j) {
      threads[j] = std::thread(computeBlock, i + j);
    }
    
    for (unsigned int j = 0; j < numThreads && i + j < numBlocks; ++j) {
      if (threads[j].joinable()) {
        threads[j].join();
      }
    }
  }
  
  
  
  // Sum up the vectors element-wise
  Eigen::VectorXd arg1_sum = Eigen::VectorXd::Zero(arg1[0].size());
  Eigen::VectorXd arg2_sum = Eigen::VectorXd::Zero(arg2[0].size());
  Eigen::VectorXd arg3_sum = Eigen::VectorXd::Zero(arg3[0].size());
  
  for (int i = 0; i < numBlocks; ++i) {
    arg1_sum += arg1[i]; 
    arg2_sum += arg2[i];
    arg3_sum += arg3[i];
  }
  
  
  
  //Perform element-wise multiplication
  Eigen::VectorXd productResult = arg2_sum.cwiseProduct(arg3_sum);
  
  //Take the square root of each element in the result
  Eigen::VectorXd sqrtResult = productResult.array().sqrt();
  
  //compute final step of the genetic correlation
  Eigen::VectorXd gencor = arg1_sum.array() / sqrtResult.array();
  
  
  return gencor;
}




// [[Rcpp::export]]

// Function to calculate heritability based on multiple simulations at once
Eigen::VectorXd calcHeritabilityManyWThreads(const Eigen::MatrixXd& sVCMat, const Eigen::MatrixXd& kinMat) {
  Eigen::MatrixXd kaVCMat = sVCMat * kinMat;
  Eigen::MatrixXd argTransposed = kaVCMat * sVCMat.transpose();
  
  Eigen::VectorXd arg = argTransposed.diagonal();
  double trKK = (kinMat * kinMat).trace();
  Eigen::VectorXd h_sq_hat(arg.size());
  
  auto computeChunk = [&](int start, int end) {
    for (int i = start; i < end; ++i) {
      h_sq_hat(i) = arg(i) / trKK;
    }
  };
  
  const unsigned int numThreads = std::thread::hardware_concurrency();
  std::vector<std::thread> threads(numThreads);
  
  int chunkSize = arg.size() / numThreads;
  int start = 0;
  
  // Create threads
  for (unsigned int i = 0; i < numThreads; ++i) {
    int end = (i == numThreads - 1) ? arg.size() : start + chunkSize;
    threads[i] = std::thread(computeChunk, start, end);
    start = end;
  }
  
  // Join threads
  for (std::thread &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }
  
  return h_sq_hat;
}