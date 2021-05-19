library(RcppArmadillo)
library(Rcpp)
sourceCpp("R/haar_cpp_openmpi.cpp")
sourceCpp("R/haar_cpp.cpp")
set.seed(12345)
library(simts)
library(microbenchmark)
size = 10e3

Xt = gen_gts(n=size, model = WN(1))
microbenchmark(wvar_cpp_openmpi(Xt, ncores = 4), wvar_cpp(Xt), times = 10)
