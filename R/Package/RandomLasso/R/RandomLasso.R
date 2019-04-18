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

RandomLasso <- function(independent, dependent, bootstraps,
                        alpha.a = 1, alpha.b = 1, verbose = FALSE) { # x & y
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

  # The cutoff value if (-cutoff < 0 < cutoff) then reduce to 0.
  cut.off <- 1 / number.of.samples

  # If argument "bootstraps" is not set, then we will set it here.
  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 80
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
    random.features <- sort(sample(features.names, number.of.samples,
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
    # random.independent.sd[[ii]] = sqrt(apply(random.independent.scale[[ii]]^2, 2, sum))
    # random.independent.sd.scale[[ii]] = scale(random.independent.scale[[ii]], FALSE, random.independent.sd[[ii]])
    # Prints time estimation.
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
  for (ii in 1:bootstraps) {
    # Sample features column numbers equal to the number of samples.
    random.features <- sort(sample(features.names, number.of.samples,
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
  # Dividing weight by number of bootstraps for final answer.
  sum.weights <- colSums(all.weights) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0
  # colnames(sum.weights) <- real.features.names
  return(sum.weights)
}
