# clean ws
rm(list = ls())

#load packages
library(wv)
library(simts)
require(magrittr)

#generate ts
set.seed(12345)
Xt = gen_gts(n=10^5, model = WN(1))

#calculate wavelet coefficient with modwt coefficient
my_haar = function(Xt){
  #import library
  require(magrittr)
  tsl = length(Xt)
  #define J
  J = log(tsl, 2) %>% floor()
  #create matrix to store elements
  haar_coeff_list = list()
  all_j =   J %>% seq()
  #define all scales
  scales = 2^all_j
  for(i_j in all_j){
    i_scale = scales[i_j]
    length_haar_transfo = tsl - 2^i_j + 1
    
    #define positive and negative yt
    coef_length = seq(length_haar_transfo)
    coef_scale_i = c()
    for(cl in coef_length){
      #define all position
      pos_scale_i = seq(cl, i_scale+cl-1)
      #define t/2 for scale j
      mid_id = length(pos_scale_i)/2
      #define left and right position
      neg_id = head(pos_scale_i, mid_id)
      pos_id = tail(pos_scale_i,length(pos_scale_i)/2)
      #calculate haar coefficient define as the weighted mean where weight equal -1, 1
      coef_scale_i = c(coef_scale_i, mean(c(-Xt[neg_id], Xt[pos_id])))
      #append to growing vector
      haar_coeff_list[[i_j]] = coef_scale_i
    }
  }
  #calculate wavelet variance
  wvariance = c(NA, length(all_j))
  for(i in seq(length(all_j))){
    # wvariance[i] = var(haar_coeff_list[[i]])*(length(haar_coeff_list[[i]])-1)/length(haar_coeff_list[[i]])
    # wvariance[i] = var(haar_coeff_list[[i]])
    wvariance[i] = t(haar_coeff_list[[i]]) %*% haar_coeff_list[[i]] / length(haar_coeff_list[[i]])
  }
  #return haar coefficients and wvariance
  return(list(haar_coeff_list, wvariance))
}


manual_computation = my_haar(Xt)

# test equality values of coefficients
for(i in seq( length(manual_computation[[1]]))) {
  myscale = 2^i
  print(paste("Coefficients equal at scale", myscale))
  print(all(round(manual_computation[[1]][[i]], 12) == round(res[[i]], 12)))
}

# test variance values
haar_coeff_list = my_haar(Xt)
cbind(round(haar_coeff_list[[2]][1:8], 10),
      round(wv::wvar(Xt)$variance, 10)
)



