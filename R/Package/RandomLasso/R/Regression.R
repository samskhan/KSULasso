Lasso <- function(independent, dependent, alpha) {
  # Alpha 1 is Lasso, Alpha 0.5 is Net Elastic, and Alpha 0 is Ridge Regression.
  cv.lasso.results <- cv.glmnet(independent, dependent, type.measure = "mse",
                                nfold = 5, alpha = alpha)
  lasso.results <- glmnet(independent, dependent, lambda = cv.lasso.results$lambda.min,
                          standardize = FALSE, alpha = alpha, intercept = FALSE)
  # Coeffiecents of random sample.
  hat.beta <- coef(lasso.results)[-1]
  return(hat.beta)
}
