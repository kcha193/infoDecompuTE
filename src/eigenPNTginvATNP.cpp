#include <RcppArmadillo.h>   

// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]   
arma::mat eigenPNTginvATNP(arma::mat& z, arma::mat& N, arma::mat& T, arma::mat& P) { 

	return arma::eig_sym(z * T * arma::trans(N) * P * N * T * z) ;
}    

