#' Performs variable selection and regularization using Random Lasso.
#'
#'
#' @param independent Matrix of treatment data.
#' @param dependent Matrix of control data.
#' @param bootstrap Number of times random sampling occures.
#' @param suppress Supresses all printing and time estimation.
#' @keywords
#' @author James Matthew Hamilton
#' @export
#' @examples
#' RandomLasso(asthma.treatment, asthma.control)

#source("src/EstimateTime.R")

RandomLasso <- function(independent, dependent, bootstrap,
                    PVAL = 0.05, C_I = 0.95, P = 0.5, S_E = 0.05,
                    REP = 10, suppress = FALSE) {
  
  sample.names <- colnames(independent)
  sample.count <- ncol(independent)
  features = nrow(independent)
  
  if (missing(bootstrap)) {
    bootstrap <- round(sample.count / features) * 80
  }
  index.and.weight = list(0,0)
  all.weights <- matrix(0, nrow = bootstrap, ncol = sample.count)
  colnames(all.weights) <- sample.names
  
  sample = list(0,0)
  for (ii in 1:bootstrap){
    sample.index <- sample(sample.names , features, replace = FALSE)
    sample[[ii]] <- independent[, which(!is.na(match(sample.names, sample.index)))]
  }
  
  print("Step 1...")
  start = as.numeric(Sys.time())
  #a = mclapply(sample, BootstrapStep1, dependent,
  #             sample.names, features, mc.cores = 2)
  a = lapply(sample, BootstrapStep1, dependent, sample.names, features)
  print(paste("Vectorized: ",as.numeric(Sys.time() - start)))
  start = as.numeric(Sys.time())
  for (jj in 1:bootstrap) {
    index.and.weight[[jj]] <- BootstrapStep1(sample[[jj]], dependent, sample.names,
                                      features)
  }
  print(paste("For-loop: ",as.numeric(Sys.time() - start)))
  for (jj in 1:bootstrap) {
    all.weights[jj, index.and.weight[[jj]][[1]]] <- index.and.weight[[jj]][[2]]
  }
   
  sum.weights <- colSums(all.weights) / bootstrap
  probability <- (1 / bootstrap) + abs(sum.weights)
  # probability <- sum.weights / sum(sum.weights)
  print("Step 2...")
  for (jj in 1:bootstrap) {
    
    index.and.weight[jj] <- BootstrapStep2(independent, dependent,  sample.names,
                                        probability,
                                        PVAL = 0.01, features)
  }
  for (jj in 1:bootstrap) {
    all.weights[jj, index.and.weight[[jj]][[1]]] <- index.and.weight[[jj]][[2]]
  }
  sum.weights <- colSums(all.weights)^2 / bootstrap
  return(sum.weights)
}

BootstrapStep1 <- function(independent, dependent, sample.names, features) {
  
  #sample.index <- sample(sample.names , features, replace = FALSE)
  #sample <- independent[, which(!is.na(match(sample.names, sample.index)))]
  
  las2 <- glmnet(independent, y, alpha = 1, family = "gaussian")
  cvlas2 <- cv.glmnet(independent, y, alpha = 1)
  hat_beta2 <- predict.glmnet(las2, type = "coefficients", s = cvlas2$lambda.1se)
  idx <- which(hat_beta2 != 0)
  
  if (length(idx) == 0) { 
    Adj_temp_b[, jj] <- 0
  } else {
    gene_id2 <- match(hat_beta2@Dimnames[[1]][-c(1:2)], sample.names)
  }
  results <- list(gene_id2, hat_beta2[-c(1,2)])
  return(results)
}

BootstrapStep2 <- function(independent, y, sample.names, probability, PVAL, features) {
  
  sample.index <- sample(sample.names, features, replace = FALSE, prob = probability)
  sample <- independent[, which(!is.na(match(sample.names, sample.index)))]
  
  las2 <- glmnet(sample, y, alpha = 1, family = "gaussian")
  cvlas2 <- cv.glmnet(sample, y, alpha = 1)
  hat_beta2 <- predict.glmnet(las2, type = "coefficients", s = cvlas2$lambda.1se)
  idx <- which(hat_beta2 != 0)
  
  if (length(idx) == 0) { 
    Adj_temp_b[, jj] <- 0
  } else {
    gene_id2 <- match(hat_beta2@Dimnames[[1]][-c(1:2)], sample.names)
  }
  results <- list(gene_id2, hat_beta2[-c(1,2)])
  return(results)
}

SampleFeatures <- function(independent, y, sample.names, probability, PVAL, features) {
  
  sample.index <- sample(sample.names, features, replace = FALSE)
  sample <- independent[, which(!is.na(match(sample.names, sample.index)))]
  
}