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

ParallelRandomLasso <- function(x, y, bootstraps, alpha = c(1, 1), verbose = TRUE, test = FALSE, cores) {

  if (test) {start = as.numeric(Sys.time())}
  x = as.matrix(x)
  y = as.matrix(y)
  features <- ncol(x)
  samples <- nrow(x)
  cut.off <- 1 / samples

  # If argument "bootstraps" is not set, then we will set it here.
  if (missing(bootstraps)) {
    bootstraps <- round(features / samples) * 40
  }

  # ------------ Step I ------------ #

  if (verbose) {
    cat("\nPart 1 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  beta.hat <- mclapply(seq_len(bootstraps), .helper.part.a, x, y,
                     features, samples,
                     bootstraps, pb, as.numeric(Sys.time()), alpha[1], verbose, mc.cores = cores)

  reduce.a <- function(x, y) {abs(x + y)}
  importance.measure <- Reduce(reduce.a, beta.hat)

  # ------------ Step II ------------ #

  if (verbose) {
    cat("\nPart 2 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  beta.hat <- mclapply(seq_len(bootstraps), .helper.part.a, x, y,
                     features, samples,
                     bootstraps, pb, as.numeric(Sys.time()), alpha[1], verbose, mc.cores = cores)

  if (verbose) {
    cat("\n [Done] \n")
    close(pb)
  }

  reduce.b <- function(x, y) {x + y / bootstraps}
  sum.weights <- Reduce(reduce.b, beta.hat)
  sum.weights[abs(sum.weights) < cut.off] <- 0
  sum.weights <- matrix(sum.weights, nrow = features, ncol = 1)
  rownames(sum.weights) <- colnames(x)
  colnames(sum.weights) <- "Coefficients"
  if (test) {
    return(c(features, samples, bootstraps, as.numeric(Sys.time()) - start))
  }
  return(sum.weights)
}
