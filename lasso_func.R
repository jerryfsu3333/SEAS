lasso_func <- function(x, y){
  cvfit <- cv.glmnet(x, y, intercept = FALSE)
  as.matrix(coef(cvfit, s = "lambda.1se"))
}