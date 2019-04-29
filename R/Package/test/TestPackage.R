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
# load("res/sim1_sig3_our.RData")
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
# >>>>>>>>>>>> 100 Features 50 Samples <<<<<<<<<<<<
x = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/qx")[,-1])
y = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/qy")[,-1])
# >>>>>>>>>>>> 1000 Features 100 Samples <<<<<<<<<<<<
x = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/x")[,-1])
y = as.matrix(read.csv("~/Dropbox/KSULasso/R/Package/test/res/y")[,-1])


detach("package:RandomLasso", unload = TRUE)
install.packages("~/Dropbox/KSULasso/R/Package/RandomLasso/", repos = NULL,
                 type = "source")
library(RandomLasso)

start <- Sys.time()
RapidRandomLasso(x, y, alpha = c(0.6, 0.6), verbose = TRUE, test = FALSE)
Sys.time() - start

start <- Sys.time()
RandomLasso(x, y, alpha = c(0.6, 0.6), verbose = TRUE, test = FALSE)
Sys.time() - start

part1 <- SteppedRandomLasso(x, y, bootstraps = 300, alpha = 1, verbose = TRUE)
part2 <- SteppedRandomLasso(x, y, importance = part1, bootstraps = 300, alpha = 1, verbose = TRUE)

# >>>>>>>>>> PERFORMANCE TESTING <<<<<<<<<<

# >>>>>> Generating Huge Test Data <<<<<<
map = list()
for (ii in 1:20) {map[[ii]] = x}
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
pb <- txtProgressBar(min = 0, max = 30, style = 3)
for (ii in 1:30) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RandomLasso(superX2[1:(ii * 20),1:(ii * 200)], superY[1:(ii * 20),], alpha = c(1, 1), verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/XY1Core", Sys.time(),".csv"))


# >>>>>> Testing Bootstraps <<<<<<
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 88, max = 200, style = 3)
for (ii in 88:200) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RandomLasso(superX2[1:50,1:100], superY[1:50,], alpha = c(1, 1), bootstraps = 5 * ii, verbose = FALSE, test = TRUE)
}
for (ii in 1:87) {
  test[[ii]] = c(0,0,0,0)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/Bootstraps(200)", Sys.time(),".csv"))

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
write.csv(x = delist,paste("../log/XYMULTICore", Sys.time(),".csv"))

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

# >>>>>> Comparing Rapid v Random Lasso Accuracy [Features: 100] <<<<<<

time1 = list()
time2 = list()
coef1 = list()
coef2 = list()
error1 = list()
error2 = list()
iterations = 30
truth <- c(3, 3, -3, 2, 2, -2, 1.5, 1.5, 1.5, -1.5)

for (ii in 1:iterations) {
  start <- Sys.time()
  coef1[[ii]] <- RandomLasso(x, y, alpha = c(0.6, 0.6), verbose = TRUE, test = FALSE)
  error1[[ii]] <- abs(coef1[[1]][1:10] - truth)
  time1[[ii]] <- Sys.time() - start
  
  start <- Sys.time()
  coef2[[ii]] <- RapidRandomLasso(x, y, alpha = c(0.6, 0.6), verbose = TRUE, test = FALSE)
  error2[[ii]] <- abs(coef2[[1]][1:10] - truth)
  # sqrt((sd(coef2[[1]][1:10])^2/10) + sd(truth)^2/10)
  time2[[ii]] <- Sys.time() - start
}
results <- matrix(NA, nrow = 10, ncol = 7)
colnames(results) <- c("Run Time Random Lasso", "Run Time Rapid Lasso", "Truth", "Average Coef Random Lasso", "Average Coef Rapid Lasso", "Difference Random Lasso", "Difference Rapid Lasso")
results[1,1] <- Reduce('+', time1) / iterations
results[1,2] <- Reduce('+', time2) / iterations
results[1:10,3] <- truth
results[1:10,4] <- (Reduce('+', coef1) / iterations)[1:10]
results[1:10,5] <- (Reduce('+', coef2) / iterations)[1:10]
results[1:10,6] <- Reduce('+', error1) / iterations
results[1:10,7] <- Reduce('+', error2) / iterations

write.csv(x = results,paste("../log//RapidAccuracyDataset1", Sys.time(),".csv"))

# >>>>>> Comparing Rapid v Random Lasso Accuracy [Features: 1000] <<<<<<

time1 = list()
time2 = list()
coef1 = list()
coef2 = list()
error1 = list()
error2 = list()
iterations = 15
truth <- c(0.5285, -2.3780, -0.8379, -1.2074, -0.7173, -0.1747, -2.4928, -1.0483, -2.7500, 2.1934, 2.9586, -1.2497, 0.8943, 4.4044, 2.2559)

for (ii in 1:iterations) {
  start <- Sys.time()
  coef1[[ii]] <- RandomLasso(x, y, alpha = c(0.6, 0.6), verbose = TRUE, test = FALSE)
  error1[[ii]] <- abs(coef1[[1]][1:15] - truth)
  time1[[ii]] <- Sys.time() - start
  
  start <- Sys.time()
  coef2[[ii]] <- RapidRandomLasso(x, y, alpha = c(0.6, 0.6), verbose = TRUE, test = FALSE)
  error2[[ii]] <- abs(coef2[[1]][1:15] - truth)
  # sqrt((sd(coef2[[1]][1:10])^2/10) + sd(truth)^2/10)
  time2[[ii]] <- Sys.time() - start
}
results <- matrix(NA, nrow = 50, ncol = 7)
colnames(results) <- c("Run Time Random Lasso", "Run Time Rapid Lasso", "Truth", "Average Coef Random Lasso", "Average Coef Rapid Lasso", "Difference Random Lasso", "Difference Rapid Lasso")
results[1,1] <- Reduce('+', time1) / iterations
results[1,2] <- Reduce('+', time2) / iterations
results[1:15,3] <- truth
results[1:50,4] <- (Reduce('+', coef1) / iterations)[1:50]
results[1:50,5] <- (Reduce('+', coef2) / iterations)[1:50]
results[1:15,6] <- Reduce('+', error1) / iterations
results[1:15,7] <- Reduce('+', error2) / iterations

write.csv(x = results,paste("../log/RapidAccuracyDataset2", Sys.time(),".csv"))

# >>>>>> Testing Rapid Random Lasso <<<<<<

test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 30, style = 3)
for (ii in 1:30) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RapidRandomLasso(superX2[1:(ii * 20),1:(ii * 200)], superY[1:(ii * 20),], alpha = c(1, 1), verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/Rapid1Core", Sys.time(),".csv"))

# >>>>>> Testing Bootstraps <<<<<<
test = list()
cat("\nTesting:\n")
pb <- txtProgressBar(min = 0, max = 200, style = 3)
for (ii in 1:200) {
  setTxtProgressBar(pb, ii)
  test[[ii]] = RapidRandomLasso(superX2[1:50,1:100], superY[1:50,], alpha = c(1, 1), bootstraps = 5 * ii, verbose = FALSE, test = TRUE)
}
delist = t(sapply(test, unlist))
write.csv(x = delist,paste("../log/RapidBootstraps", Sys.time(),".csv"))








