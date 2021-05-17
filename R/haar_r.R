# rm(list=ls())
# library(simts)

#calculate wavelet coefficient with modwt coefficient
wvar_r = function(Xt){
  tsl = length(Xt)
  #define J
  J = log(tsl, 2) %>% floor()
  #create list to store elements
  haar_coeff_list = list()
  all_j =   J %>% seq()
  #define all scales
  scales = 2^all_j
  for(i_j in all_j){
    i_scale = scales[i_j]
    length_haar_transfo = tsl - 2^i_j + 1
    
    #define positive and negative yt
    coef_length = seq(length_haar_transfo)
    coef_scale_i = vector(length =length_haar_transfo, mode = "numeric")
    for(cl in coef_length){
      #define all position
      pos_scale_i = seq(cl, i_scale+cl-1)
      #define t/2 for scale j
      mid_id = length(pos_scale_i)/2
      #define left and right position
      pos_id = pos_scale_i[1:mid_id] 
      neg_id =  tail(pos_scale_i,length(pos_scale_i)/2)
      #calculate haar coefficient define as the weighted mean where weight equal -1, 1
      xt_neg = Xt[neg_id] * -1
      xt_pos = Xt[pos_id] 
      xt_weighted = c(xt_neg, xt_pos)
      coef_scale_i[cl] = mean(xt_weighted)
    }
    #append to growing vector
    haar_coeff_list[[i_j]] = coef_scale_i
  }
  #calculate wavelet variance
  wvariance = vector(mode = "numeric", length = length(all_j))
  for(i in seq(length(all_j))){
    # wvariance[i] = var(haar_coeff_list[[i]])*(length(haar_coeff_list[[i]])-1)/length(haar_coeff_list[[i]])
    # wvariance[i] = var(haar_coeff_list[[i]])
    haar_coef = haar_coeff_list[[i]]
    wvariance[i] = t(haar_coef) %*% haar_coef / length(haar_coef)
  }
  #return haar coefficients and wvariance
  return(wvariance)
}



# Xt = gen_gts(n=10^3, model = WN(1))
# wv::wvar(Xt)
# wvar_r(Xt)
