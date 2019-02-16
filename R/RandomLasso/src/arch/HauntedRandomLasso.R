#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' @param independent Matrix of treatment data.
#' @param dependent Matrix of control data.
#' @param bootstrap Number of times random sampling occures.
#' @param suppress Supresses all printing and time estimation.
#' @keywords
#' @author James Matthew Hamilton
#' @export
#' @examples
#' RandomLasso(asthma.treatment, asthma.control)

#source("src/EstimateTime.R")

RandomLasso <- function(independent, dependent, bootstrap, suppress = FALSE) {

  #independent <- Normalize(independent)
  #dependent <- Normalize(dependent)
  sample.names <- colnames(independent)
  sample.count <- ncol(independent)
  features <- nrow(independent)
  lambda <- 10^seq(10, -2, length = 100)
  sample = list(0,0)
  sample.index = list(0,0)
  index.and.beta = list(0,0)
  
  if (missing(bootstrap)) {
    bootstrap <- round(sample.count / features) * 20
  }
  all.weights <- matrix(0, nrow = bootstrap, ncol = sample.count)
  colnames(all.weights) <- sample.names
  
  for (ii in 1:bootstrap) {
    sample.index[[ii]] <- sample(1:sample.count, features, replace = FALSE)
    sample.index[[ii]] <- sort(sample.index[[ii]])
    sample[[ii]] <- independent[, which(!is.na(match((1:sample.count),
                                                     sample.index[[ii]])))]
  }
  print("Step 1......")
  #start = as.numeric(Sys.time())
  #a = mclapply(sample, BootstrapStep1, dependent,
  #              sample.names, features, mc.cores = 2)
  #a = lapply(sample, BootstrapStep1, dependent, sample.names, features, lambda)
  #print(paste("Vectorized: ",as.numeric(Sys.time() - start)))
  #start = as.numeric(Sys.time())
  for (jj in 1:bootstrap) {
    index.and.beta[[jj]] <- BootstrapStep1(sample[[jj]], dependent, lambda)
  }
  #print(paste("For-loop: ",as.numeric(Sys.time() - start)))
  for (jj in 1:bootstrap) {
    all.weights[jj, sample.index[[jj]]] <- index.and.beta[[jj]]
  }
  sum.weights <- colSums(all.weights) / bootstrap
  probability <- (1 / bootstrap) + abs(sum.weights)
  # probability <- sum.weights / sum(sum.weights)
  for (ii in 1:bootstrap) {
    sample.index[[ii]] <- sample(1:sample.count, features, replace = FALSE,
                                 prob = probability)
    sample.index[[ii]] <- sort(sample.index[[ii]])
    sample[[ii]] <- independent[, which(!is.na(match((1:sample.count),
                                                     sample.index[[ii]])))]
  }
  print("Step 2...")
  for (jj in 1:bootstrap) {
    index.and.beta[jj] <- BootstrapStep2(sample[[jj]], dependent, lambda)
  }
  for (jj in 1:bootstrap) {
    all.weights[jj, sample.index[[jj]]] <- index.and.beta[[jj]]
  }
  sum.weights <- colSums(all.weights) / bootstrap
  return(sum.weights)
}

BootstrapStep1 <- function(independent, dependent, lambda) {
  lasso.results <- glmnet(independent, dependent, alpha = 1, lambda = lambda)
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients",
                             s = cv.lasso.results$lambda.min)
  index.and.beta <- hat.beta[-c(1)]
  return(index.and.beta)
}

BootstrapStep2 <- function(independent, dependent, lambda) {
  lasso.results <- glmnet(independent, dependent, alpha = 1, lambda = lambda)
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients",
                             s = cv.lasso.results$lambda.min)
  index.and.beta <- hat.beta[-c(1)]
  return(index.and.beta)
}

Normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
SampleFeatures <- function(independent, dependent, sample.names, probability, PVAL, features) {
  
  sample.index <- sample(sample.names, features, replace = FALSE)
  sample <- independent[, which(!is.na(match(sample.names, sample.index)))]
  
}