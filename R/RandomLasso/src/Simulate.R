## Simulation DataSet Explanation  ##


#Data type : Rdata
#Dataset name :  sim1_sig3_our.RData / sim2_sig3_our.RData  / sim3_sig3_our.RData / sim4_sig3_our.RData 
#i : iteration (from 1 to 10)

#training set variables :  x, y  
#validataion set variables : x.val, y.val
#testing set variables : x.test, y.test
#covariance matrix : cov_x

#R program fro reading data

setwd("~/KSULasso/R/RandomLasso/")
load("res/sim1_sig3_our.RData")
#sim_data2 = load("res/sim2_sig3_our.RData")
#sim_data3 = load("res/sim3_sig3_our.RData")
#sim_data4 = load("res/sim4_sig3_our.RData")
#sim_data = sim_data2

n_iter <- 1    # iteration

for (i in 1:n_iter)   # i is index of replication
{
  # load simulation data
  y <- sim_data[[i]][2]
  y <- unlist(y)
  y <- matrix(y, ncol = 1)
  
  x <- sim_data[[i]][3]
  x <- unlist(x[1])
  x <- matrix(x, ncol = n_feature, byrow = FALSE)
  
  y.test <- sim_data[[i]][4]
  y.test <- unlist(y.test)
  y.test <- matrix(y.test, ncol = 1)
  
  x.test <- sim_data[[i]][5]
  x.test <- unlist(x.test[1])
  x.test <- matrix(x.test, ncol = n_feature, byrow = FALSE)
  
  y.val  <- sim_data[[i]][6]
  y.val <- unlist(y.val)
  y.val <- matrix(y.val, ncol = 1)
  
  x.val  <- sim_data[[i]][7]
  x.val <- unlist(x.val[1])
  x.val <- matrix(x.val, ncol = n_feature, byrow = FALSE)
  
  cov_x <- sim_data[[i]][8]   # same sigmax
  cov_x <- unlist(cov_x[1])
  cov_x <- matrix(cov_x, ncol = n_feature, byrow = FALSE)
  
  # x.test, y.test, = loadData(sim_data)
  
  n_feature <- ncol(x)   
  colnames(x)      <- paste('P', seq(1:n_feature), sep='')
  colnames(x.val)  <- paste('P', seq(1:n_feature), sep='')
  colnames(x.test) <- paste('P', seq(1:n_feature), sep='')
  
}

setwd("~/KSULasso/R/RandomLasso/")
source("src/RandomLasso.R")
lasso.coef <- RandomLasso(x, y)