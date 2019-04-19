.helper.part.a <- function(ii, x, y, number.of.features,
                           number.of.samples, bootstraps, pb,
                           start.time, alpha, verbose, importance) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

  # Sample features column numbers equal to the number of samples.

  random.features <- sample(number.of.features, number.of.samples,
                            replace = FALSE)
  # Mix up the rows.
  random.samples <- sample(number.of.samples, replace = TRUE)
  # Subset the columns and rows from the independent data.
  random.x <- x[random.samples, random.features]
  # Subset the rows from the dependent data.
  random.y <- y[random.samples, ]

  # Centering the dependent variable.
  random.y.mean <- mean(random.y)
  random.y.scale <- random.y - random.y.mean

  # Centering the independent variable.
  random.x.mean <- apply(random.x, 2, mean)
  random.x.scale <- scale(random.x, random.x.mean, FALSE)
  standard.deviation <- sqrt(apply(random.x.scale^2, 2, sum))
  random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)
  # Filling in the empty rows one by one with results from lasso.
  beta.hat[random.features] <- Lasso(random.x.scale,random.y.scale,
                                     alpha)
  return(beta.hat)
}

.helper.part.b <- function(ii, x, y, number.of.features,
                           number.of.samples, bootstraps, pb,
                           start.time, alpha, verbose, importance) {
  beta.hat <- numeric(number.of.features)
  if (verbose) {.helper.time.remaining(pb, start.time, ii, bootstraps)}

  # Sample features column numbers equal to the number of samples.
  random.features <- sample(number.of.features, number.of.samples,
                            replace = FALSE, prob = importance)
  # Mix up the rows.
  random.samples <- sample(number.of.samples, replace = TRUE)
  # Subset the columns and rows from the independent data.
  random.x <- x[random.samples, random.features]
  # Subset the rows from the dependent data.
  random.y <- y[random.samples, ]

  # Centering the dependent variable.
  random.y.mean <- mean(random.y)
  random.y.scale <- random.y - random.y.mean

  # Centering the independent variable.
  random.x.mean <- apply(random.x, 2, mean)
  random.x.scale <- scale(random.x, random.x.mean, FALSE)
  standard.deviation <- sqrt(apply(random.x.scale^2, 2, sum))
  random.x.scale <- scale(random.x.scale, FALSE, standard.deviation)

  # Filling in the empty rows one by one with results from lasso.
  beta.hat[random.features] <- Lasso(random.x.scale,random.y.scale,
                                     alpha) / standard.deviation
  return(beta.hat)
}
