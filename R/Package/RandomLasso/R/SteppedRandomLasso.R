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
#' SteppedRandomLasso(independent, dependent)

SteppedRandomLasso <- function(x, y, importance, bootstraps, alpha = 1, verbose = TRUE) {

  x = as.matrix(x)
  y = as.matrix(y)
  number.of.features <- ncol(x)
  number.of.samples <- nrow(x)

  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 80
  }

  if (missing(importance)) {
    if (verbose) {
      cat("Running Part 1 Only:\n")
      pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
    }
    beta.hat <- lapply(seq_len(bootstraps), .helper.part.a, x, y,
                     number.of.features, number.of.samples,
                     bootstraps, pb, as.numeric(Sys.time()), alpha, verbose)
  } else {
    importance = abs(as.matrix(importance))
    if (verbose) {
      cat("Running Part 2 Only:\n")
      pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
    }
    beta.hat <- lapply(seq_len(bootstraps), .helper.part.b, x, y,
                       number.of.features, number.of.samples,
                       bootstraps, pb, as.numeric(Sys.time()), alpha, verbose,
                       importance)
  }

  if (verbose) {
    cat("\n [Done]")
    close(pb)
  }

  sum.weights <- Reduce('+', beta.hat) / bootstraps
  sum.weights <- matrix(sum.weights, nrow = number.of.features, ncol = 1)
  rownames(sum.weights) <- colnames(x)
  colnames(sum.weights) <- "Coefficients"
  return(sum.weights)
}
