#' Performs variable selection and regularization using Random Lasso and Spark.
#'
#'
#' @param independent Matrix of intependent data.
#' @param dependent Matrix of dependent data.
#' @param bootstraps Number of times random sampling occurs.
#' @param suppress Supresses all printing and time estimation.
#' @param cutoff A cut off value for reducing insignificant values to zero.
#' @keywords
#' @author James Matthew Hamilton
#' @export
#' @examples
#' RandomLasso(independent, dependent)


RandomLassoSp <- function(independent, dependent, bootstraps, suppress = FALSE,
                        cutoff = TRUE) { # x & y
  library(glmnet)
  library(sparklyr)
  source("src/EstimateTime.R")
  config <- spark_config()
  config$spark.executor.cores <- 2
  config$spark.executor.memory <- "4G"
  sc <- spark_connect(master = "local")
  # We save the real column names, since we are replacing with numeric.
  real.features.names <- colnames(independent)
  number.of.features <- ncol(independent)
  number.of.samples <- nrow(independent)
  # Replacing row and column names with numeric values.
  # This arrays of numeric column names will come in handy later.
  features.names <- (1:number.of.features) # Array of column names.
  sample.names <- (1:number.of.samples) # Array of row names.
  colnames(independent) <- features.names
  row.names(independent) <- sample.names
  row.names(dependent) <- sample.names
  # Declaring rest of variables into memory.
  random.features <- list()
  random.independent <- list()
  random.dependent <- list()
  # The cutoff value if (-cutoff < 0 < cutoff) then reduce to 0.
  cut.off <- 1 / number.of.samples
  
  # If argument "bootstraps" is not set, then we will set it here.
  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 100
  }
  # Declairing an empty matrix.
  all.weights <- matrix(0, nrow = bootstraps, ncol = number.of.features)
  colnames(all.weights) <- features.names
  
  ##########
  # Step I #
  ##########
  # For-loop that bootstraps samples. This can be easily vectorized and
  # parallelized, but we are waitng to finish implementation
  for (ii in 1:bootstraps) {
    # Sample features column numbers equal to the number of samples.
    random.features[[ii]] <- sort(sample(features.names, number.of.samples,
                                         replace = FALSE))
    # Mix up the rows.
    random.samples <- sample(number.of.samples, replace = TRUE)
    # Subset the columns and rows from the independent data.
    random.independent[[ii]] <- independent[random.samples, random.features[[ii]]]
    # Subset the rows from the dependent data.
    random.dependent[[ii]] <- dependent[random.samples, ]
  }
  # Start timer, will add feature to supress printouts later.
  # This is in function argument "supress".
  start = Sys.time()
  print("Step 1..")
  # This is anouther easy area to make vectorized and parallelized,
  # i.e. spark_apply()
  # a = lapply(sample, BootstrapStep1, dependent, sample.names, features)
  sdf_len(sc, random.independent) %>%
    spark_apply(Coefficients, random.dependent)
  print(a)
  break
  for (ii in 1:bootstraps) {
    # Prints time estimation.
    if (ii %% 20 == 0) {PrintEstimateTime(start, ii, bootstraps)}
    # Filling in the empty rows one by one with results from lasso.
    all.weights[ii, random.features[[ii]]]  <- Coefficients(random.independent[[ii]],
                                                            random.dependent[[ii]])
  }
  # Getting the sum of ever column.
  importance.measure <- (1 / 1000000000) + abs(colSums(all.weights))
  
  ###########
  # Step II #
  ###########
  for (ii in 1:bootstraps) {
    # Using the sum of every column to add probabilty weight for next sampling.
    random.features[[ii]] <- sort(sample(features.names, number.of.samples,
                                         replace = FALSE, prob = importance.measure))
    # Rest is same as before...
    random.samples <- sample(number.of.samples, replace = TRUE)
    random.independent[[ii]] <- independent[random.samples, random.features[[ii]]]
    random.dependent[[ii]] <- dependent[random.samples, ]
  }
  start = Sys.time()
  print("Step 2...")
  for (ii in 1:bootstraps) {
    if (ii %% 20 == 0) {PrintEstimateTime(start, ii, bootstraps)}
    all.weights[ii, random.features[[ii]]] <- Coefficients(random.independent[[ii]],
                                                           random.dependent[[ii]])
  }
  # Dividing weight by number of bootstraps for final answer.
  sum.weights <- colSums(all.weights) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0
  #colnames(sum.weights) <- real.features.names
  print(sum.weights)
  return(sum.weights)
}

##############################
# Lasso Coeffiecent Function #
##############################

Coefficients <- function(independent, dependent) {
  
  # Alpha 1 is Lasso, Alpha 0.5 is Net Elastic, and Alpha 0 is Ridge Regression.
  lasso.results <- glmnet(independent, dependent, alpha = 1, family = "gaussian")
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  # Coeffiecents of random sample.
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients",
                             s = cv.lasso.results$lambda.1se)
  #print(hat.beta)
  return(hat.beta[-1])
}