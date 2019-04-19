#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' @param x Matrix of intependent data.
#' @param y Matrix of dependent data.
#' @param bootstraps Number of random bootstraps taken.
#' @param alpha A cut off value for reducing insignificant values to zero.
#' @param verbose Supresses all printing and time estimating functions.
#' @keywords
#' @author James Matthew Hamilton
#' @export
#' @examples
#' RandomLasso(x, y)
#' RandomLasso(x, y, alpha.a = 0.9, alpha.b = 0.9, verbose = FALSE, bootstraps = 300)
#'

RandomLasso <- function(x, y, bootstraps, alpha = c(1, 1), verbose = TRUE) {

  x = as.matrix(x)
  y = as.matrix(y)
  number.of.features <- ncol(x)
  number.of.samples <- nrow(x)
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

  beta.hat <- lapply(seq_len(bootstraps), .helper.part.a, x, y,
         number.of.features, number.of.samples,
         bootstraps, pb, as.numeric(Sys.time()), alpha[1], verbose)

  importance.measure <- abs(Reduce('+', beta.hat))

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
  return(sum.weights)
}
