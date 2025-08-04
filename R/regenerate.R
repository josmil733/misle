#' @export regenerate

regenerate <- function(r,y,X,theta,beta,A){
      m <- A %*% y
      n <- length(m)
      predictor <- X%*%theta
      trans_pr <- exp(predictor+beta*m)/(exp(predictor+beta*m)+exp(-predictor-beta*m))
      return(list(y=2*rbinom(n, 1, trans_pr)-1))
}
