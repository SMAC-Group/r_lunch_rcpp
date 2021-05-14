#include <Rcpp.h>
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

// [[Rcpp::export]]
NumericVector create_empty_vector(int n) {
  NumericVector v (n);
  return v;
}
// [[Rcpp::export]]
int get_length_vec(NumericVector x) {
  int length_vec = x.length();
  return length_vec;
}

// [[Rcpp::export]]
void print_i_cpp(int i) {
  for(int j = 1; j<=i; j ++){
    Rcout << j  << "\n";
  }
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
a = create_empty_vector(10)
get_length_vec(a)
print_i_cpp(3)
*/
