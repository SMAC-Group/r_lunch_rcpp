



vector(mode = "numeric", length = 5)


library(microbenchmark)
big_list = list()
# fill list
for(i in seq(1000)){
  big_list[[i]] = rnorm(10000)
}
f1_vapply = function(big_list){
  out = vapply(big_list, FUN = mean, FUN.VALUE = numeric(1))
  return(out)
}
f2_for_loop_append = function(big_list){
  out = c()
  for(i in seq(length(big_list))){
    res = mean(big_list[[i]])
    out = c(out, res)
  }
  return(out)
}

f3_for_loop= function(big_list){
  out = vector(mode = "numeric", length = length(big_list))
  for(i in seq(length(big_list))){
    out[i]  = mean(big_list[[i]])
  }
  return(out)
}



microbenchmark(f1_vapply(big_list),
               f2_for_loop_append(big_list), 
               f3_for_loop(big_list))



microbenchmark(vapply(big_list, FUN = mean, FUN.VALUE = numeric(1)),
               for(i in seq(length(big_list))){
                 mean(big_list[[i]])
               })


big_mat = matrix(rnorm(100000), ncol = 1000, nrow=100)
f1_for_loop = function(big_mat){
  out = vector(mode = "numeric", length = dim(big_mat)[2])
  for(i_col in seq(dim(big_mat)[2])){
    out[i_col] = sum(big_mat[,i_col] )
  }
  return(out)
}
f1_for_loop(big_mat)

f2_apply = function(big_mat){
  out = apply(big_mat, MARGIN = 2, FUN = function(x){sum(x)}) 
  return(out)
}
f2_apply(big_mat)

microbenchmark(f1_for_loop(big_mat),
               f2_apply(big_mat))



library(DescTools)
profvis::profvis({
  n = 150
  prob_p <- 0.1 # true proportion
  nbsim <- 1000 # nbr of simulations
  cover <- vector(mode = "numeric", length = nbsim) 
  set.seed(123)
  for (i in 1:nbsim) {
    X <- rbinom(n, 1, prob_p) # generate random variable from bernouilli
    p_hat = mean(X)
    var_p_hat = sqrt((p_hat*(1-p_hat))/n)
    ci = c( p_hat - qnorm(.975)*var_p_hat,  p_hat + qnorm(.975)*var_p_hat)
    cover[i] <- ifelse(ci[1] <= prob_p & ci[2] >= prob_p, 1, 0) # save coverage
  }
  mean(cover)
})



n = 10000L
x = seq(from = 0, to = pi, length.out = n)

f = function(x){
  
}


baseline = median(microbenchmark(
  
  y <- f(x)
  
)$time)


apply_time = median(microbenchmark(
  
  y2 <- sapply(x, f)
  
)$time)


for_time = median(microbenchmark({
  
  y3 <- vector(mode = "list", length = n)
  for(i in 1:n){
    y3[[i]] = f(x[i])
  }
  y3 <- as.numeric(y3)
  
})$time)


apply_time / baseline

for_time / baseline

for_time / apply_time
