lasso_func <- function(x, y, ...){
  # cvfit <- cv.glmnet(x, y, intercept = FALSE)
  cvfit <- cv.glmnet(x, y, ...)
  as.matrix(coef(cvfit, s = "lambda.1se"))
}