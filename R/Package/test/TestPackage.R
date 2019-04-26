## Simulation DataSet Explanation  ##


#Data type : Rdata
#Dataset name :  sim1_sig3_our.RData / sim2_sig3_our.RData  / sim3_sig3_our.RData / sim4_sig3_our.RData
#i : iteration (from 1 to 10)

#training set variables :  x, y
#validataion set variables : x.val, y.val
#testing set variables : x.test, y.test
#covariance matrix : cov_x

options(scipen = 999)

setwd("~/Dropbox/KSULasso/R/Package/test")
load("res/sim1_sig3_our.RData")
# load("res/sim2_sig3_our.RData")
# load("res/sim3_sig3_our.RData")
# load("res/sim4_sig3_our.RData")

n_iter <- 1

for (i in 1:n_iter) {
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

  cov_x <- sim_data[[i]][8]
  cov_x <- unlist(cov_x[1])
  cov_x <- matrix(cov_x, ncol = n_feature, byrow = FALSE)

  n_feature <- ncol(x)
  colnames(x)      <- paste('P', seq(1:n_feature), sep = '')
  colnames(x.val)  <- paste('P', seq(1:n_feature), sep = '')
  colnames(x.test) <- paste('P', seq(1:n_feature), sep = '')
}

rm(RandomLasso)
setwd("~/Dropbox/KSULasso/R/Package/test")
x = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/qx")[,-1])
y = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/qy")[,-1])
# x = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/x")[,-1])
# y = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/y")[,-1])


detach("package:RandomLasso", unload = TRUE)
install.packages("~/Dropbox/KSULasso/R/Package/RandomLasso/", repos = NULL,
                 type = "source")
library(RandomLasso)

RandomLasso(x, y, alpha = c(0.5, 0.5), verbose = TRUE, test = FALSE)

part1 <- SteppedRandomLasso(x, y, bootstraps = 300, alpha = 1, verbose = TRUE)
part2 <- SteppedRandomLasso(x, y, importance = part1, bootstraps = 300, alpha = 1, verbose = TRUE)

# >>>>>>>>>> PERFORMANCE TESTING <<<<<<<<<<

# >>>>>> Generating Huge Test Data <<<<<<
map = list()
for (ii in 1:10) {map[[ii]] = x}
superX = Reduce(cbind, map)
dim(superX)

mapx = list()
mapy = list()
for (ii in 1:20) {mapx[[ii]] = superX}
for (ii in 1:20) {mapy[[ii]] = y}
superY = Reduce(rbind, mapy)
superX2 = Reduce(rbind, mapx)
dim(superY)
dim(superX2)

# >>>>>> Testing X and Y <<<<<<
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 200, style = 3)
for (ii in 1:200) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RandomLasso(superX2[1:(ii * 10),1:(ii * 100)], superY[1:(ii * 10),], alpha = c(1, 1), bootstraps = 5 , verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/XY1Core", Sys.time(),".csv"))


# >>>>>> Testing Bootstraps <<<<<<
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 200, style = 3)
for (ii in 1:100) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RandomLasso(superX2[1:50,1:100], superY[1:50,], alpha = c(1, 1), bootstraps = 5 * ii, verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/Bootstraps", Sys.time(),".csv"))

# >>>>>> Parallel Testing X and Y <<<<<<
library(parallel)
detectCores()
cores <- 2
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 200, style = 3)
for (ii in 1:200) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = ParallelRandomLasso(superX2[1:(ii * 10),1:(ii * 100)], superY[1:(ii * 10),], alpha = c(1, 1), bootstraps = 5 , verbose = FALSE, test = TRUE, cores = cores)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/XY2Core", Sys.time(),".csv"))

# >>>>>> Testing Just Features <<<<<<
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 200, style = 3)
for (ii in 1:100) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RandomLasso(superX2[1:100,1:(200 * ii)], superY[1:100,], alpha = c(1, 1), bootstraps = 10, verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/Features", Sys.time(),".csv"))

# >>>>>> Testing Just Samples <<<<<<
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 200, style = 3)
for (ii in 1:100) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RandomLasso(superX2[1:(10 * ii),1:2000], superY[1:(10 * ii),], alpha = c(1, 1), bootstraps = 10, verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/Samples", Sys.time(),".csv"))
