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

RapidRandomLasso <- function(x, y, bootstraps, alpha = c(1, 1), verbose = TRUE, test = FALSE) {

  if (test) {start = as.numeric(Sys.time())}
  number.of.features <- ncol(x)
  number.of.samples <- nrow(x)
  cut.off <- 1 / number.of.samples

  # If argument "bootstraps" is not set, then we will set it here.
  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 40
  }

  if(alpha[1] == alpha[2]) {
    single.lambda <- cv.glmnet(x, y, type.measure = "mse", nfold = 10, alpha = alpha[1])$lambda.min
    lambda <- c(single.lambda, single.lambda)
    print(lambda)
  } else{
    lambda <- c(cv.glmnet(x, y, type.measure = "mse", nfold = 10, alpha = alpha[1])$lambda.min,
                cv.glmnet(x, y, type.measure = "mse", nfold = 10, alpha = alpha[2])$lambda.min)
  }

  if (verbose) {
    cat("\nPart 1 of 2:\n")
    pb <- txtProgressBar(min = 0, max = bootstraps, style = 3)
  }

  .part1 <- function(ii, x, y, start.time) {
    if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

    beta.hat <- numeric(number.of.features)
    random.features <- sample(number.of.features, number.of.samples, replace = FALSE)
    random.samples <- sample(number.of.samples, replace = TRUE)

    random.x <- x[random.samples, random.features]
    random.y <- y[random.samples, ]

    random.y.mean <- mean(random.y)
    random.y.scale <- random.y - random.y.mean

    random.x.mean <- apply(random.x, 2, mean)
    random.x.scale <- scale(random.x, random.x.mean, FALSE)
    standard.deviation <- sqrt(apply(random.x.scale^2, 2, sum))
    random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)

    beta.hat[random.features] <- RapidLasso(random.x.scale, random.y.scale, alpha[1], lambda[1])
    return(beta.hat)
  }

  beta.hat <- lapply(seq_len(bootstraps), .part1, x, y, as.numeric(Sys.time()))
  importance.measure <- abs(Reduce('+', beta.hat))
  print(importance.measure)

  # ------------ Step II ------------ #

  .part2 <- function(ii, x, y, start.time) {
    if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

    beta.hat <- numeric(number.of.features)

    random.features <- sample(number.of.features, number.of.samples, replace = FALSE, prob = importance.measure)
    random.samples <- sample(number.of.samples, replace = TRUE)

    random.x <- x[random.samples, random.features]
    random.y <- y[random.samples, ]

    random.y.mean <- mean(random.y)
    random.y.scale <- random.y - random.y.mean

    random.x.mean <- apply(random.x, 2, mean)
    random.x.scale <- scale(random.x, random.x.mean, FALSE)
    standard.deviation <- sqrt(apply(random.x.scale^2, 2, sum))
    random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)

    beta.hat[random.features] <- RapidLasso(random.x.scale, random.y.scale, alpha[2], lambda[2]) / standard.deviation

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
  sum.weights <- matrix(sum.weights, nrow = number.of.features, ncol = 1)
  rownames(sum.weights) <- colnames(x)
  colnames(sum.weights) <- "Coefficients"
  if (test) {
    return(c(number.of.features, number.of.samples, bootstraps, as.numeric(Sys.time()) - start))
  }
  return(sum.weights)
}
