library(RcppArmadillo)
library(Rcpp)
sourceCpp("R/haar_cpp_openmpi.cpp")
library(simts)
library(microbenchmark)


size = 10^1
set.seed(12345)
Xt = gen_gts(n=size, model = WN(1))
res_benchmark_size_10 = print(microbenchmark(
  wvar_cpp_openmpi(Xt, ncores = 4),
  times = 5
), unit = "s")

save(res_benchmark_size_10, file = "data/res_benchmark_size_10.rda")


size = 10^2
set.seed(12345)
Xt = gen_gts(n=size, model = WN(1))
res_benchmark_size_100 = print(microbenchmark(
  wvar_cpp_openmpi(Xt, ncores = 4),
  times = 5
), unit= "s")

save(res_benchmark_size_100, file = "data/res_benchmark_size_100.rda")


size = 10^3
set.seed(12345)
Xt = gen_gts(n=size, model = WN(1))
res_benchmark_size_1000 = print(microbenchmark(
  wvar_cpp_openmpi(Xt, ncores = 4),
  times = 5
), unit ="s")

save(res_benchmark_size_1000, file = "data/res_benchmark_size_1000.rda")


size = 10^4
set.seed(12345)
Xt = gen_gts(n=size, model = WN(1))
res_benchmark_size_10000 = print(microbenchmark(
  wvar_cpp_openmpi(Xt, ncores = 4),
  times = 5
), unit = "s")

save(res_benchmark_size_10000, file = "data/res_benchmark_size_10000.rda")


size = 10^5
set.seed(12345)
Xt = gen_gts(n=size, model = WN(1))
res_benchmark_size_100000 = print(microbenchmark(
  wvar_cpp_openmpi(Xt, ncores = 4),
  times = 5
), unit ="s")
beepr::beep()

save(res_benchmark_size_100000, file = "data/res_benchmark_size_100000.rda")


bench_open_mpi = rbind(
  res_benchmark_size_10,
  res_benchmark_size_100,
  res_benchmark_size_1000,
  res_benchmark_size_10000,
  res_benchmark_size_100000
)
save(bench_open_mpi, file = "data/bench_open_mpi.rda")

bench_open_mpi
mat_mean_perf_2


load("data/mat_mean_perf.rda")

mat_mean_perf_2 = mat_mean_perf /1e6
par(mar=c(5,7,1,1))
plot(x = 10^seq(5), y = mat_mean_perf_2[,2], 
     log="x", type = "b", xlab = "Sample size", ylim =c(0, max(mat_mean_perf_2)),
     cex.lab = 1.2, las = 1, 
     cex.main = 1.5,cex.axis = 1.1,
     col = "darkblue", pch = 16, yaxt ="n",ylab="")
mtext(side = 2, text = "Mean execution time (seconds)", line = 5.6, cex= 1.1)
axis(side=2, las = 2, at = )
grid(col="grey80", lty=1)
polygon(x = c(10^seq(5), rev(10^seq(5))), y = c(mat_mean_perf_2[,3], rev(mat_mean_perf_2[,5])), col = "#FF8C0019" , border = NA)
polygon(x = c(10^seq(5), rev(10^seq(5))), y = c(mat_mean_perf_2[,4], rev(mat_mean_perf_2[,6])), col = "#00008B19", border = NA)
lines(x = 10^seq(5), y = mat_mean_perf_2[,1], type ="b",  col = "darkorange", pch = 16)
lines(x = 10^seq(5), y = mat_mean_perf_2[,2], type ="b",  col = "darkblue", pch = 16)
lines(x = 10^seq(5), y = bench_open_mpi$mean)


legend("topleft",col = c("darkorange", "darkblue"), legend= c("R", "CPP"), lwd=1, pch = 16, bty ="n", cex = 1.5)
