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

RandomLasso <- function(independent, dependent, bootstraps,
                        alpha.a = 1, alpha.b = 1, verbose = TRUE) {

  number.of.features <- ncol(independent)
  number.of.samples <- nrow(independent)
  cut.off <- 1 / number.of.samples

  # If argument "bootstraps" is not set, then we will set it here.
  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 80
  }

  # ------------ Step I ------------ #

  if (verbose) {
  cat("Part 1 of 2:\n")
  pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  ii = 0
  beta.hat = lapply(seq_len(bootstraps), .helper.part.a, independent, dependent,
         number.of.features, number.of.samples,
         bootstraps, pb, as.numeric(Sys.time()), alpha.a, verbose)

  importance.measure <- abs(Reduce('+', beta.hat))

# ------------ Step II ------------ #

  if (verbose) {
    cat("\nPart 2 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  ii = 0
  beta.hat = lapply(seq_len(bootstraps), .helper.part.b, independent, dependent,
             number.of.features, number.of.samples,
             bootstraps, pb, as.numeric(Sys.time()), alpha.b, verbose,
             importance.measure)

  if (verbose) {
    cat("\n [Done]")
    close(pb)
  }

  sum.weights <- Reduce('+', beta.hat) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0.0
  sum.weights <- matrix(sum.weights, nrow = number.of.features, ncol = 1)
  rownames(sum.weights) <- colnames(independent)
  colnames(sum.weights) <- "Coefficients"
  return(sum.weights)
}

.helper.part.a <- function(ii, independent, dependent, number.of.features,
                                    number.of.samples,bootstraps, pb,
                                    start.time, alpha, verbose) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

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
  beta.hat[random.features] <- Lasso(random.independent.scale,random.dependent.scale,
                                     alpha) / standard.deviation
  return(beta.hat)
}

.helper.part.b <- function(ii, independent, dependent, number.of.features,
                           number.of.samples,bootstraps, pb,
                           start.time, alpha, verbose, importance) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

  # Sample features column numbers equal to the number of samples.
  random.features <- sort(sample(number.of.features, number.of.samples,
                                 replace = FALSE, prob = importance))
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
  beta.hat[random.features] <- Lasso(random.independent.scale,random.dependent.scale,
                                     alpha) / standard.deviation
  return(beta.hat)
}

.helper.time.remaining <- function(pb, start.time, current.increment,
                                   end.increment) {
  setTxtProgressBar(pb, current.increment)
  passed <- as.numeric(Sys.time()) - start.time
  remaining <- (passed / (current.increment / end.increment)) - passed
  hr <- floor(remaining / 3600)
  min <- floor(remaining / 60) - (hr * 60)
  sec <- floor(remaining) - (hr * 3600) - (min * 60)

  cat("\r", paste("[",hr, ":", min, ":", sec, "] |", sep = ""))
}
