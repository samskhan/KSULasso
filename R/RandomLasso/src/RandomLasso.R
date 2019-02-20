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

RandomLasso <- function(independent, dependent, bootstrap, suppress = FALSE,
                        cutoff = TRUE) {
  real.feature.names <- colnames(independent)
  features <- ncol(independent)
  samples <- nrow(independent)
  feature.index <- (1:features)
  index.feature <- (1:samples)
  colnames(independent) <- feature.index
  row.names(independent) <- index.feature
  row.names(dependent) <- index.feature
  cut.off <- 1 / samples
  index.and.beta = list(0,0)
  sample.independent = list(0,0)
  sample.dependent = list(0,0)
  if (missing(bootstrap)) {
    bootstrap <- round(features / samples) * 40
  }
  all.weights <- matrix(0, nrow = bootstrap, ncol = features)
  colnames(all.weights) <- feature.index


  for (ii in 1:bootstrap) {
    index.feature <- sample(feature.index, samples, replace = FALSE)
    index.sample <- sample(samples, samples, replace = TRUE)
    sample.independent[[ii]] <- independent[index.sample, which(!is.na(match(feature.index, index.feature)))]
    sample.dependent[[ii]] <- dependent[index.sample, ]
  }
  print("Step 1..")
  for (ii in 1:bootstrap) {
    index.and.beta[[ii]] <- Lasso(sample.independent[[ii]], sample.dependent[[ii]], feature.index)
  }
  for (ii in 1:bootstrap) {
    all.weights[ii, index.and.beta[[ii]][[1]]] <- index.and.beta[[ii]][[2]]
  }

  importance.measure <- (1 / 1000000000) + abs(colSums(all.weights))
  for (ii in 1:bootstrap) {
    index.feature <- sample(feature.index, samples, replace = FALSE, prob = importance.measure)
    index.sample <- sample(samples, samples, replace = TRUE)
    sample.independent[[ii]] <- independent[index.sample, which(!is.na(match(feature.index, index.feature)))]
    sample.dependent[[ii]] <- dependent[index.sample, ]
  }
  print("Step 2...")
  for (ii in 1:bootstrap) {
    index.and.beta[[ii]] <- Lasso(sample.independent[[ii]], sample.dependent[[ii]], feature.index)
  }
  for (ii in 1:bootstrap) {
    all.weights[ii, index.and.beta[[ii]][[1]]] <- index.and.beta[[ii]][[2]]
  }
  sum.weights <- colSums(all.weights) / bootstrap
  sum.weights[sum.weights < cut.off] <- 0
  print(sum.weights)
  return(sum.weights)
}

Lasso <- function(independent, dependent, sample.names) {
  
  lasso.results <- glmnet(independent, dependent, alpha = 1, family = "gaussian")
  cv.lasso.results <- cv.glmnet(independent, dependent, alpha = 1)
  hat.beta <- predict.glmnet(lasso.results, type = "coefficients",
                             s = cv.lasso.results$lambda.1se)
  print(hat.beta)
  gene.index <- match(hat.beta@Dimnames[[1]][-1], sample.names)
  print(gene.index)
  index.and.beta <- list(gene.index, hat.beta[-1])
  return(index.and.beta)
}