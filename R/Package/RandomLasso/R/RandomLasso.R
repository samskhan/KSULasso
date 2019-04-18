#' Performs variable selection and regularization using Random Lasso.
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
#'

.helper.time.remaining <- function(pb, start.time, current.increment,
                                         end.increment){
  setTxtProgressBar(pb, current.increment)
  passed <- as.numeric(Sys.time()) - start.time
  remaining <- (passed / (current.increment / end.increment)) - passed
  hr <- floor(remaining / 3600)
  min <- floor(remaining / 60) - (hr * 60)
  sec <- floor(remaining) - (hr * 3600) - (min * 60)

  cat("\r", paste("[",hr, ":", min, ":", sec, "] |", sep = ""))
  flush.console()

}

RandomLasso <- function(independent, dependent, bootstraps,
                        alpha.a = 1, alpha.b = 1, verbose = TRUE) {

  number.of.features <- ncol(independent)
  number.of.samples <- nrow(independent)
  cut.off <- 1 / number.of.samples

  # If argument "bootstraps" is not set, then we will set it here.
  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 80
  }
  # Declairing an empty matrix.
  all.weights <- matrix(0, nrow = bootstraps, ncol = number.of.features)
  colnames(all.weights) <- colnames(independent)

  ##########
  # Step I #
  ##########
  # For-loop that bootstraps samples. This can be easily vectorized and
  # parallelized, but we are waitng to finish implementation
  if (verbose) {
  cat("Part 1 of 2:\n")
  start <- as.numeric(Sys.time())
  pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }
  for (ii in 1:bootstraps) {
    if (verbose) {.helper.time.remaining(pb, start, ii, bootstraps)}
    # Sample features column numbers equal to the number of samples.
    random.features <- sort(sample(number.of.features, number.of.samples,
                                         replace = FALSE))
    # Mix up the rows.
    random.samples <- sample(number.of.samples, replace = TRUE)
    # Subset the columns and rows from the independent data.
    random.independent <- independent[random.samples, random.features]
    # Subset the rows from the dependent data.
    random.dependent <- dependent[random.samples, ]

    # Centering the dependent variable.
    random.dependent.mean <- mean(random.dependent)
    random.dependent.scale <- random.dependent - random.dependent.mean

    # Centering the independent variable.
    random.independent.mean <- apply(random.independent, 2, mean)
    random.independent.scale <- scale(random.independent, random.independent.mean, FALSE)
    standard.deviation <- sqrt(apply(random.independent.scale^2, 2, sum))
    random.independent.scale <- scale(random.independent.scale, FALSE, standard.deviation)

    # Obtaining the standard deviation.
    # Filling in the empty rows one by one with results from lasso.
    all.weights[ii, random.features] <- Lasso(random.independent.scale,
                                                    random.dependent.scale,
                                                    alpha.a)
  }
  # Getting the sum of ever column.
  importance.measure <- abs(colSums(all.weights))

  ###########
  # Step II #
  ###########
  if (verbose) {
    cat("\nPart 2 of 2:\n")
    start <- as.numeric(Sys.time())
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }
  for (ii in 1:bootstraps) {
    if (verbose) {.helper.time.remaining(pb, start, ii, bootstraps)}
    # Sample features column numbers equal to the number of samples.
    random.features <- sort(sample(number.of.features, number.of.samples,
                                   replace = FALSE, prob = importance.measure))
    # Mix up the rows.
    random.samples <- sample(number.of.samples, replace = TRUE)
    # Subset the columns and rows from the independent data.
    random.independent <- independent[random.samples, random.features]
    # Subset the rows from the dependent data.
    random.dependent <- dependent[random.samples, ]

    # Centering the dependent variable.
    random.dependent.mean <- mean(random.dependent)
    random.dependent.scale <- random.dependent - random.dependent.mean

    # Centering the independent variable.
    random.independent.mean <- apply(random.independent, 2, mean)
    random.independent.scale <- scale(random.independent, random.independent.mean, FALSE)
    standard.deviation <- sqrt(apply(random.independent.scale^2, 2, sum))
    random.independent.scale <- scale(random.independent.scale, FALSE, standard.deviation)
    # Filling in the empty rows one by one with results from lasso.
    all.weights[ii, random.features] <- Lasso(random.independent.scale,
                                              random.dependent.scale,
                                              alpha.b) / standard.deviation
  }
  close(pb)
  # Dividing weight by number of bootstraps for final answer.
  sum.weights <- colSums(all.weights) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0
  # colnames(sum.weights) <- real.features.names
  return(sum.weights)
}
