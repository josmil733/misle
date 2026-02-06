#' @export step1.mnh

step1.mnh <- function(grid.lambda,
                         d,
                         eps,
                         y.mvc,
                         y = if (sample.split) y.is.train else y.is,
                         X = if (sample.split) X.is.train else X.is,
                         beta,
                         A.m,
                         mvc.c,
                         tau = 0.8,
                         theta){

require(magrittr)

N <- length(y)

    lam <- seq(grid.lambda$from, grid.lambda$to, length.out = grid.lambda$length.out)
    coef_seq <- matrix(0, nrow = d, ncol = length(lam))
    # proximal gradient method
    for (s in 1:length(lam)) {
      gam_ok <- rep(0, d)
      gam_kk <- rep(100, d)
      t_k <- 1
      k <- 0
      while (sum(abs(gam_ok-gam_kk))>eps) {
        gam_k <- gam_ok
        grad <- lossgrad(
          theta = gam_k,
          beta = beta,
          A.m = A.m,
          mvc.c = mvc.c,
          y.mvc = y.mvc,
          y.is = y,
          X.is = X
        )
        proxi <- prox(
          gradL = grad,
          theta = gam_k,
          lam = lam[s],
          t = t_k
        )
        G <- (gam_k-proxi)/(t_k)
        ## line search
        if (loss(
          theta = gam_k - t_k * G,
          beta = beta,
          A.m = A.m,
          mvc.c = mvc.c,
          y.mvc = y.mvc,
          y.is = y,
          X.is = X)
          > loss(
            theta = gam_k,
            beta = beta,
            A.m = A.m,
            mvc.c = mvc.c,
            y.mvc = y.mvc,
            y.is = y,
            X.is = X)
            -t_k*sum(grad*(G))+t_k*sum(G**2)/2){
          t_k <- t_k*tau
        }else{
          gam_ok <- proxi
          gam_kk <- gam_k
          k <- k+1
        }
      }
      coef_seq[, s] <- gam_ok
    }

    # select model with BIC: 2*L_N+log(N)*no_nonzero/N
    BIC <- numeric(length(lam))
    log_like <- numeric(length(lam))
    for (k in 1:length(lam)) {
        log_like[k] <- loss(
          theta = coef_seq[,k],
          beta = beta,
          A.m = A.m,
          mvc.c = mvc.c,
          y.mvc = y.mvc,
          y.is = y,
          X.is = X
        )
      BIC[k] <- 2*log_like[k] + log(N)*length(which(as.vector(coef_seq[, k]) != 0))/N
    }

    min_indx <- which.min(BIC)

    theta_mincv <- coef_seq[, min_indx]

    ## Define final items for output

    results = matrix(NA, nrow = length(lam), ncol = 6)

    indices.on <- which(theta != 0) #only used for comparison at the end. This value does not inform the method.

    return(list(
      coef = theta_mincv,
      best.lambda = lam[min_indx],
      # results = data.frame(lambda=lambda, ind.1 = coef_seq[indices.on[1], ], ind.2 = coef_seq[indices.on[2], ], count.on=colSums(coef_seq != 0), BIC = BIC)
      results = matrix(
        c(
          lam %>% round(-log(lam[1], 10)+1),
          coef_seq[indices.on[1], ],
          coef_seq[indices.on[2], ],
          colSums(coef_seq != 0),
          BIC,
          apply(coef_seq, 2, subtract, theta) %>% abs() %>% colSums()
        )
      )
    ))

}

