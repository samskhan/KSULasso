#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' @param independent Matrix of treatment data.
#' @param dependent Matrix of control data.
#' @param bootstraps Number of times random sampling occures.
#' @param suppress Supresses all printing and time estimation.
#' @keywords
#' @author James Matthew Hamilton
#' @export
#' @examples
#' RandomLasso(independent, dependent)

library(sparklyr)
source("src/EstimateTime.R")


RandomLasso <- function(independent, dependent, bootstraps, suppress = FALSE,
                        cutoff = TRUE) {
  real.features.names <- colnames(independent)
  number.of.features <- ncol(independent)
  number.of.samples <- nrow(independent)
  features.names <- (1:number.of.features)
  column.names <- (1:number.of.samples)
  colnames(independent) <- features.names
  row.names(independent) <- column.names
  row.names(dependent) <- column.names
  random.features <- list()
  random.independent <- list()
  random.dependent <- list()
  beta.hat <- list()
  cut.off <- 1 / number.of.samples

  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 100
  }
  all.weights <- matrix(0, nrow = bootstraps, ncol = number.of.features)
  colnames(all.weights) <- features.names

  for (ii in 1:bootstraps) {
    random.features[[ii]] <- sort(sample(features.names, number.of.samples,
                                         replace = FALSE))
    random.samples <- sample(number.of.samples, replace = TRUE)
    random.independent[[ii]] <- independent[random.samples, random.features[[ii]]]
    random.dependent[[ii]] <- dependent[random.samples, ]
  }
  start = Sys.time()
  print("Step 1..")
  for (ii in 1:bootstraps) {
    if (ii %% 20 == 0) {PrintEstimateTime(start, ii, bootstraps)}
    beta.hat[[ii]] <- Lasso(random.independent[[ii]], random.dependent[[ii]])
  }
  for (ii in 1:bootstraps) {
    all.weights[ii, random.features[[ii]]] <- beta.hat[[ii]]
  }

  importance.measure <- (1 / 1000000000) + abs(colSums(all.weights))
  
  for (ii in 1:bootstraps) {
    random.features[[ii]] <- sort(sample(features.names, number.of.samples,
                                          replace = FALSE, prob = importance.measure))
    random.samples <- sample(number.of.samples, replace = TRUE)
    random.independent[[ii]] <- independent[random.samples, random.features[[ii]]]
    random.dependent[[ii]] <- dependent[random.samples, ]
  }
  start = Sys.time()
  print("Step 2...")
  for (ii in 1:bootstraps) {
    if (ii %% 20 == 0) {PrintEstimateTime(start, ii, bootstraps)}
    beta.hat[[ii]] <- Lasso(random.independent[[ii]], random.dependent[[ii]])
  }
  for (ii in 1:bootstraps) {
    all.weights[ii, random.features[[ii]]] <- beta.hat[[ii]]
  }
  sum.weights <- colSums(all.weights) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0
  #colnames(sum.weights) <- real.features.names
  print(sum.weights)
  return(sum.weights)
}

Lasso <- function(independent, dependent) {
  
  lasso.results <- glmnet(independent, dependent, alpha = 1, family = "gaussian")
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients",
                             s = cv.lasso.results$lambda.1se)
  #print(hat.beta)
  return(hat.beta[-1])
}