#' @export regenerate

regenerate <- function(r, y) {
      m <- A %*% y
      X <- matrix(runif(n * d, min = -k, max = k), nrow = n, ncol = d)
      predictor <- X %*% theta
      trans_pr <- exp(predictor + beta * m) /
            (exp(predictor + beta * m) + exp(-predictor - beta * m))
      return(list(y = 2 * rbinom(n, 1, trans_pr) - 1, X = X))
}
