Lasso <- function(x, y, alpha) {
  # Alpha 1 is Lasso, Alpha 0.5 is Net Elastic, and Alpha 0 is Ridge Regression.
  cv.lasso.results <- cv.glmnet(x, y, type.measure = "mse",
                                nfold = 5, alpha = alpha)
  lasso.results <- glmnet(x, y, lambda = cv.lasso.results$lambda.min,
                          standardize = FALSE, alpha = alpha, intercept = FALSE)
  # Coeffiecents of random sample.
  hat.beta <- coef(lasso.results)[-1]
  return(hat.beta)
}
