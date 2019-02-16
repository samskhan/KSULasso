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

RandomLasso <- function(independent, dependent, bootstrap, suppress = FALSE) {
  
  sample.names <- colnames(independent)
  features <- ncol(independent)
  samples = nrow(independent)
  index.and.beta = list(0,0)
  sample = list(0,0)
  if (missing(bootstrap)) {
    bootstrap <- round(features / samples) * 80
  }
  all.weights <- matrix(0, nrow = bootstrap, ncol = features)
  colnames(all.weights) <- sample.names

  for (ii in 1:bootstrap) {
    sample.index <- sample(sample.names , samples, replace = FALSE)
    sample[[ii]] <- independent[, which(!is.na(match(sample.names, sample.index)))]
  }
  print("Step 1..")
  for (jj in 1:bootstrap) {
    index.and.beta[[jj]] <- BootstrapStep1(sample[[jj]], dependent, sample.names)
  }
  for (jj in 1:bootstrap) {
    all.weights[jj, index.and.beta[[jj]][[1]]] <- index.and.beta[[jj]][[2]]
  }
  sum.weights <- colSums(all.weights) / bootstrap
  probability <- (1 / bootstrap) + abs(sum.weights)
  # probability <- sum.weights / sum(sum.weights)
  for (ii in 1:bootstrap) {
    sample.index <- sample(sample.names, samples, replace = FALSE, prob = probability)
    sample[[ii]] <- independent[, which(!is.na(match(sample.names, sample.index)))]
  }
  print("Step 2...")
  for (jj in 1:bootstrap) {
    index.and.beta[jj] <- BootstrapStep2(sample[[jj]], dependent, sample.names)
  }
  for (jj in 1:bootstrap) {
    all.weights[jj, index.and.beta[[jj]][[1]]] <- index.and.beta[[jj]][[2]]
  }
  sum.weights <- colSums(all.weights) / bootstrap
  return(sum.weights)
}

BootstrapStep1 <- function(independent, dependent, sample.names) {
  
  lasso.results <- glmnet(independent, dependent, alpha = 1, family = "gaussian")
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients", s = cv.lasso.results$lambda.1se)
  gene.index <- match(hat.beta@Dimnames[[1]][-1], sample.names)
  index.and.beta <- list(gene.index, hat.beta[-1])
  return(index.and.beta)
}

BootstrapStep2 <- function(independent, dependent, sample.names) {
  
  lasso.results <- glmnet(independent, dependent, alpha = 1, family = "gaussian")
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients", s = cv.lasso.results$lambda.1se)
  gene.index <- match(hat.beta@Dimnames[[1]][-1], sample.names)
  index.and.beta <- list(gene.index, hat.beta[-1])
  return(index.and.beta)
}

SampleFeatures <- function(independent, dependent, sample.names, probability, samples) {
  
  sample.index <- sample(sample.names, samples, replace = FALSE)
  sample <- independent[, which(!is.na(match(sample.names, sample.index)))]
}