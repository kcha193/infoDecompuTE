#include <RcppArmadillo.h>   

// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]   
arma::mat PNTginvATNP(arma::mat& z, arma::mat& N, arma::mat& T, arma::mat& P) { 

	return z * N * T * P * T * arma::trans(N) *arma::trans(z) ;
}    

