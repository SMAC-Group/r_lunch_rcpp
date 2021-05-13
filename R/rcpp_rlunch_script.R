#########################################
# Working with cpp in R
#########################################
# clean ws
rm(list = ls())

#load packages
library(wv)
library(simts)
require(magrittr)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)

################# measuring performance in R


# loading functions with source cpp
sourceCpp("R/sum_cpp.cpp")

# loading function with inline cpp



##################################
# Basic exemple
##################################





##################################
# Example Wavelet Variance
##################################

# load R implementation
sourceCpp(file = "R/haar_cpp.cpp")
# load cpp implementation
source(file = "R/haar_r.R")

# test equality
#generate ts
set.seed(12345)
Xt = gen_gts(n=10^3, model = WN(1))

# compute wvar
wvar_r_implementation = wvar_r(Xt)
wvar_cpp_implementation = wvar_cpp(Xt)
wvar_package_implementation = wv::wvar(Xt)
all.equal(as.vector(wvar_cpp_implementation), wvar_r_implementation)
all.equal(wvar_r_implementation[1:8], wvar_package_implementation$variance)



# benchmark for different signal size
lst_benchmark = list()
for(i in seq(4)){
  size = 10^i
  set.seed(12345)
  Xt = gen_gts(n=size, model = WN(1))
  res_benchmark_size_i = microbenchmark(
    wvar_r(Xt),
    wvar_cpp(Xt),
    times = 10
  )
  lst_benchmark[[i]] = print(res_benchmark_size_i, unit ="us")
  print(size)
}
lst_benchmark

mat_mean_perf = matrix(ncol = 6, nrow= 4)
for(i in seq(4)){
  mat_mean_perf[i, ] = cbind(lst_benchmark[[i]][,4], lst_benchmark[[i]][,2],lst_benchmark[[i]][,7])
}

getwd()
save(mat_mean_perf, file = "data/mat_mean_perf.rda")

t_col <- function(color, percent = 90, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  t.col
}
t_col("darkblue")
t_col("darkorange", percent = .8)

par(mar=c(5,5,1,1))
plot(x = 10^seq(4), y = mat_mean_perf[,2], 
     log="x", type = "b", ylim = c(1, 10e6), xlab = "Sample size",
     ylab = "Mean execution time (microseconds)", cex.lab = 1.5, cex.main = 1.5,cex.axis = 1.3,
     col = "darkblue", pch = 16)
grid(col="grey80", lty=1)
polygon(x = c(10^seq(4), rev(10^seq(4))), y = c(mat_mean_perf[,3], rev(mat_mean_perf[,5])), col = "#FF8C0019" , border = NA)
polygon(x = c(10^seq(4), rev(10^seq(4))), y = c(mat_mean_perf[,4], rev(mat_mean_perf[,6])), col = "#00008B19", border = NA)
lines(x = 10^seq(4), y = mat_mean_perf[,1], type ="b",  col = "darkorange", pch = 16)
lines(x = 10^seq(4), y = mat_mean_perf[,2], type ="b",  col = "darkblue", pch = 16)
legend("topleft",col = c("darkorange", "darkblue"), legend= c("R", "CPP"), lwd=1, pch = 16, bty ="n", cex = 1.5)




