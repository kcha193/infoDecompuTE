#include <RcppArmadillo.h>   

// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]   
Rcpp::List eigenRcpp(arma::mat& M) { 
	arma::vec eigval;
	arma::mat eigvec;

	eig_sym(eigval, eigvec, M);     

	return Rcpp::List::create(Rcpp::Named("value") = eigval,
                          Rcpp::Named("vector") = eigvec);
}    

