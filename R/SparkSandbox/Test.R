library(sparklyr)
library(dplyr)

spark_install()
conf <- spark_config()
conf$sparklyr.cores.local <- 2
conf$sparklyr.shell.driver-memory <- "8G"
conf$spark.memory.fraction <- 0.9

sc <- spark_connect(master = "local", config = conf)

setwd("~/Dropbox/KSULasso/R/SparkSandbox/")
x = read.csv("x", row.names = 1)
y = read.csv("y", row.names = 1)
colnames(y) <- "dependent"
all <- cbind(y, x)
alpha = 1
verbose = TRUE
number.of.features <- ncol(x)
number.of.samples <- nrow(x)
bootstraps <- round(number.of.features / number.of.samples) * 80

all <- copy_to(sc, all, overwrite = TRUE)
empty <- copy_to(sc, all, overwrite = TRUE)

# map (random features, ___)
#
# Coef
# emit()
#
# Reduce()
partition = list()

for (ii in 1:bootstraps) {
partition[[ii]] = all %>%
  select(c(1, sample(number.of.features, number.of.samples, replace = FALSE))) %>%
  sdf_sample(replacement = TRUE) %>%
  ml_linear_regression(response = dependent ~ .,
                      fit_intercept = TRUE,
                      lastic_net_param = 1)
}

# >>>>>>>>>> Spark Test <<<<<<<<<<
if (verbose) {
  cat("Spark Part:\n")
  pb <- txtProgressBar(min = 0, max = 10, style = 3)
}
spark.time <- Sys.time()
spark.beta.hat <- lapply(seq_len(10), SparkBootstrap, all,
                  number.of.features, number.of.samples,
                  10, pb, as.numeric(Sys.time()), alpha, verbose)
spark.time <- Sys.time() - spark.time
print(spark.time)



# >>>>>>>>>> Normal R Test <<<<<<<<<<
if (verbose) {
  cat("Normal R Part:\n")
  pb <- txtProgressBar(min = 0, max = 10, style = 3)
}
r.time <- Sys.time()
r.beta.hat <- lapply(seq_len(10), RBootstrap, x, y,
                   number.of.features, number.of.samples,
                   10, pb, as.numeric(Sys.time()), alpha, verbose)
r.time <- Sys.time() - r.time
print(r.time)


# >>>>>>>>>> FUNCTIONS <<<<<<<<<<
SparkBootstrap <- function(ii, all, number.of.features,
                           number.of.samples, bootstraps, pb,
                           start.time, alpha, verbose, importance) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

  random.features <- sample(number.of.features, number.of.samples,
                            replace = FALSE)
  random.all <- all %>%
    select(c(1,random.features)) %>%
    sdf_sample(replacement = FALSE)

  fit <- random.all %>% ml_linear_regression(response = dependent ~ .,
                                             fit_intercept = TRUE,
                                             lastic_net_param = 1)
  return(fit)
}

RBootstrap <- function(ii, x, y, number.of.features,
                           number.of.samples, bootstraps, pb,
                           start.time, alpha, verbose) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

  # Sample features column numbers equal to the number of samples.
  random.features <- sample(number.of.features, number.of.samples,
                            replace = FALSE)
  # Mix up the rows.
  random.samples <- sample(number.of.samples, replace = TRUE)
  # Subset the columns and rows from the independent data.
  random.x <- x[random.samples, random.features]
  # Subset the rows from the dependent data.
  random.y <- y[random.samples, ]

  # Centering the dependent variable.
  random.y.mean <- mean(random.y)
  random.y.scale <- random.y - random.y.mean

  # Centering the independent variable.
  random.x.mean <- apply(random.x, 2, mean)
  random.x.scale <- scale(random.x, random.x.mean, FALSE)
  standard.deviation <- sqrt(apply(random.x.scale^2, 2, sum))
  random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)

  # Filling in the empty rows one by one with results from lasso.
  beta.hat[random.features] <- Lasso(random.x.scale, random.y.scale, alpha) / standard.deviation

  cv.lasso.results <- cv.glmnet(random.x.scale, random.y.scale, type.measure = "mse",
                                nfold = 5, alpha = alpha)
  lasso.results <- glmnet(random.x.scale, random.y.scale, lambda = cv.lasso.results$lambda.min,
                          standardize = FALSE, alpha = alpha, intercept = FALSE)
  # Coeffiecents of random sample.
  hat.beta <- coef(lasso.results)[-1]
  return(hat.beta)
}

