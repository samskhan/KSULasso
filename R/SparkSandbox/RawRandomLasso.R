# >>>>>>>>>>>> MAIN <<<<<<<<<<<<

  install.packages("glmnet")
  install.packages("dplyr")
  install.packages("sparklyr")

  library(glmnet)
  library(sparklyr)
  library(dplyr)
  # sessionInfo() # to view package versions

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
  alpha = c(0.5, 0.5)
  verbose = TRUE
  number.of.features <- ncol(x)
  number.of.samples <- nrow(x)
  bootstraps <- round(number.of.features / number.of.samples) * 80
  cut.off <- 1 / number.of.samples
  empty <- matrix(NA, ncol = number.of.features, nrow = bootstraps)

  # x_sp <- copy_to(sc, x, overwrite = TRUE)
  # y_sp <- copy_to(sc, y, overwrite = TRUE)
  all <- copy_to(sc, all, overwrite = TRUE)
  empty <- copy_to(sc, empty, overwrite = TRUE)


# ------------ Playground ------------

  x_sp %>% head()
  x_sp %>% colnames()

  partitions <- x_sp %>%
    select(c(1,3,9))

  partitions %>% head()

  avgs <- summarize_all(x_sp, funs(mean)) %>% as.data.frame()
  exprs <- as.list(paste(colnames(x_sp),"-", avgs))
  x_sp %>%
    spark_dataframe() %>%
    invoke("selectExpr", exprs) %>%
    invoke("toDF", as.list(colnames(x_sp))) %>%
    invoke("registerTempTable", "centered")
  xx <- tbl(sc, "centered")
  xx %>% head()

  random.features <- sample(number.of.features,
                            number.of.samples,
                            replace = FALSE)
  random.all <- all %>%
    select(c(1,random.features))

  random.all %>% ml_linear_regression(response = dependent ~ .,
                                             fit_intercept = TRUE,
                                             lastic_net_param = 1)

  iris_tbl <- sdf_copy_to(sc, iris, overwrite = TRUE)
  iris_tbl %>% spark_apply(function(e){cv.glmnet(as.matrix(x), as.matrix(y),
                                                 type.measure = "mse",
                                                 nfold = 5, alpha = 1)$lambda})

  all %>% nrow2()

  # ------------ Playground END ------------

  # ------------ Step I ------------ #

  if (verbose) {
    cat("Part 1 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  beta.hat3 <- lapply(seq_len(2), .helper.part.a, all,
                     number.of.features, number.of.samples,
                     bootstraps, pb, as.numeric(Sys.time()), alpha[1], verbose)

  importance.measure <- abs(Reduce('+', beta.hat3))

  # ------------ Step II ------------ #

  if (verbose) {
    cat("\nPart 2 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  beta.hat <- lapply(seq_len(bootstraps), .helper.part.b, x, y,
                     number.of.features, number.of.samples,
                     bootstraps, pb, as.numeric(Sys.time()), alpha[2], verbose,
                     importance.measure)

  if (verbose) {
    cat("\n [Done]")
    close(pb)
  }

  sum.weights <- Reduce('+', beta.hat) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0
  sum.weights <- matrix(sum.weights, nrow = number.of.features, ncol = 1)
  rownames(sum.weights) <- colnames(x)
  colnames(sum.weights) <- "Coefficients"
  print(sum.weights)

# >>>>>>>>>>>> END OF MAIN <<<<<<<<<<<<


# >>>>>>>>>>>> FUNCTION PART A <<<<<<<<<<<<
# [Description] Randomly samples features from x and runs Random Lasso.
# [RETURNS] Coefficients from Random Lasso.
.helper.part.a <- function(ii, all, number.of.features,
                           number.of.samples, bootstraps, pb,
                           start.time, alpha, verbose, importance) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

  # Sample features column numbers equal to the number of samples.
  random.features <- sample(number.of.features, number.of.samples,
                            replace = FALSE)
  #
  random.all <- all %>%
    select(c(1,random.features))

  fit <- random.all %>% ml_linear_regression(response = dependent ~ .,
                                             fit_intercept = TRUE,
                                             lastic_net_param = 1)

  return(fit)
}

# >>>>>>>>>>>> FUNCTION PART B <<<<<<<<<<<<
# [Description] Randomly samples features WITH WEIGHT from x and runs Random Lasso.
# [RETURNS] Coefficients from Random Lasso.
.helper.part.b <- function(ii, x, y, number.of.features,
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

  return(beta.hat)
}


# >>>>>>>>>>>> FUNCTION LASSO <<<<<<<<<<<<
# [Description] Used in .helper.part.a and .helper.part.b
# [RETURNS] Coefficients from Random Lasso.
Lasso <- function(independent, dependent, alpha) {
  # Alpha 1 is Lasso, Alpha 0.5 is Net Elastic, and Alpha 0 is Ridge Regression.
  cv.lasso.results <- cv.glmnet(independent, dependent, type.measure = "mse",
                                nfold = 5, alpha = alpha)
  lasso.results <- glmnet(independent, dependent, lambda = cv.lasso.results$lambda.min,
                          standardize = FALSE, alpha = alpha, intercept = FALSE)
  # Coeffiecents of random sample.
  hat.beta <- coef(lasso.results)[-1]
  return(hat.beta)
}


# >>>>>>>>>>>> FUNCTION TIME HELPER <<<<<<<<<<<<
# [Description] Controls printing of progress bar and time remaining.
# [Returns] Nothing, prints to consol.
.helper.time.remaining <- function(pb, start.time, current.increment,
                                   end.increment) {
  setTxtProgressBar(pb, current.increment)
  passed <- as.numeric(Sys.time()) - start.time
  remaining <- (passed / (current.increment / end.increment)) - passed
  hr <- floor(remaining / 3600)
  min <- floor(remaining / 60) - (hr * 60)
  sec <- floor(remaining) - (hr * 3600) - (min * 60)
  if (hr < 10 && hr > 0) {hr = paste("0", hr, sep = "")}
  if (min < 10 && min > 0) {min = paste("0", min, sep = "")}
  if (sec < 10 && sec > 0) {sec = paste("0", sec, sep = "")}

  cat("\r", paste("[",hr, ":", min, ":", sec, "] |", sep = ""))
}
