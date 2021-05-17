#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// 


// [[Rcpp::export]]
arma::colvec wvar_cpp(arma::vec Xt) {
  int tsl = Xt.n_elem;
  // define J
  int J = floor(log(tsl) / log(2));
  List haar_coef_list = List::create();
  IntegerVector all_j = seq_len( J );
  for(int i_j = 1; i_j <= J; i_j++) {
    // Define scale tau
    double i_scale = pow(2,i_j);
    int length_haar_transfo = tsl - i_scale + 1;
    //  create empty vector to fill
    NumericVector coef_scale_i (length_haar_transfo);
    for(int cl = 1; cl <= length_haar_transfo; cl++) {
      arma::vec pos_scale_i = arma::linspace(cl-1, i_scale +cl -2,  i_scale );
      //  define negative, positive and mid id
      int mid_id = pos_scale_i.n_elem/2 ;
      arma::vec neg_id = pos_scale_i.tail(mid_id);
      arma::vec pos_id = pos_scale_i.head(mid_id);
      // Convert to position to arma uvec
      arma::uvec pos_id_2 = arma::conv_to<arma::uvec>::from(pos_id);
      arma::uvec neg_id_2 = arma::conv_to<arma::uvec>::from(neg_id);
      // extract from vector
      arma::vec xt_neg =  Xt.elem(neg_id_2) * -1;
      arma::vec xt_pos =  Xt.elem(pos_id_2) ;
      //  compute coefficient
      arma::vec xt_weighted = join_cols(xt_neg, xt_pos);
      coef_scale_i(cl-1) = mean(xt_weighted);
    }
    // append to rcpp list
    haar_coef_list.push_back(coef_scale_i);
    }
  //  create empty vector for wavelet variance
  arma::colvec wvariance (J);
  // populate wvariance with wavelet variance 
  // computed on coefficient for each scale j
  for(int i = 0; i < J; i++) {
    arma::vec haar_coef = haar_coef_list[i];
    wvariance.row(i) = haar_coef.t() * haar_coef / haar_coef.n_elem;
  }
  return wvariance;
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R


*/
