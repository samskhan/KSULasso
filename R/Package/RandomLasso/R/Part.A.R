Part.A <- function(independent, dependent, bootstraps, alpha) {

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

  if (missing(bootstraps)) {
    bootstraps <- round(number.of.features / number.of.samples) * 80
  }

  all.weights <- matrix(0, nrow = bootstraps, ncol = number.of.features)
  colnames(all.weights) <- features.names

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
    # Filling in the empty rows one by one with results from lasso.
    all.weights[ii, random.features] <- Lasso(random.independent.scale,
                                              random.dependent.scale,
                                              alpha)
  }
  return(all.weights)
}
