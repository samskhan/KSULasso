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

RandomLasso <- function(x, y, bootstraps, alpha = c(1, 1), verbose = TRUE, test = FALSE) {

  x = as.matrix(x)
  y = as.matrix(y)
  if (test) {start = as.numeric(Sys.time())}
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

  .part1 <- function(ii, x, y, start_time) {
    if (verbose) {.helper.time.remaining(pb, start_time, ii, bootstraps)}

    beta.hat <- replicate(features, 0)

    random.features <- sample(features, samples, replace = FALSE)
    random.samples <- sample(samples, replace = TRUE)

    random.x <- x[random.samples, random.features]
    random.y <- y[random.samples, ]

    random.y.mean <- mean(random.y)
    random.y.scale <- random.y - random.y.mean

    random.x.mean <- apply(random.x, 2, mean)
    random.x.scale <- scale(random.x, random.x.mean, FALSE)
    standard.deviation <- sqrt(apply(random.x.scale ^ 2, 2, sum))
    random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)

    beta.hat[random.features] <- Lasso(random.x.scale, random.y.scale,
                                       alpha[1]) / standard.deviation
    return(beta.hat)
  }

  beta.hat <- lapply(seq_len(bootstraps), .part1, x, y, as.numeric(Sys.time()))
  importance.measure <- abs(Reduce('+', beta.hat))

# ------------ Step II ------------ #

  .part2 <- function(ii, x, y, start_time) {
    if (verbose) {
      .helper.time.remaining(pb, start_time, ii, bootstraps)
    }

    beta.hat <- replicate(features, 0)

    random.features <- sample(features, samples, replace = FALSE, prob = importance.measure)
    random.samples <- sample(samples, replace = TRUE)

    random.x <- x[random.samples, random.features]
    random.y <- y[random.samples, ]

    random.y.mean <- mean(random.y)
    random.y.scale <- random.y - random.y.mean

    random.x.mean <- apply(random.x, 2, mean)
    random.x.scale <- scale(random.x, random.x.mean, FALSE)
    standard.deviation <- sqrt(apply(random.x.scale ^ 2, 2, sum))
    random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)

    beta.hat[random.features] <- Lasso(random.x.scale, random.y.scale,
                                       alpha[2]) / standard.deviation
    return(beta.hat)
  }

  if (verbose) {
    cat("\nPart 2 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  beta.hat <- lapply(seq_len(bootstraps), .part2, x, y, as.numeric(Sys.time()))

  if (verbose) {
    cat("\n [Done] \n")
    close(pb)
  }

  sum.weights <- Reduce('+', beta.hat) / bootstraps
  sum.weights[abs(sum.weights) < cut.off] <- 0
  sum.weights <- matrix(sum.weights, nrow = features, ncol = 1)
  rownames(sum.weights) <- colnames(x)
  colnames(sum.weights) <- "Coefficients"
  if (test) {
    return(c(features, samples, bootstraps, as.numeric(Sys.time()) - start))
  }
  return(sum.weights)
}
