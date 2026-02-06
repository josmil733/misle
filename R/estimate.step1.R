#' @export estimate.step1

estimate.step1 <- function(
   grid.lambda = list(from = 0.01, to = 0.1, length.out = 20),
   d,
   eps,
   y.mvc,
   y = if (sample.split) y.is.train else y.is,
   X = if (sample.split) X.is.train else X.is,
   beta,
   A.m,
   mvc.c,
   tau = 0.8,
   theta,
   max.iter=300,
   alpha.bic=alpha.bic
) {
   # attach(grid.lambda)
   n = nrow(X)
   lambda <- seq(
      from = grid.lambda$from,
      to = grid.lambda$to,
      length.out = grid.lambda$length.out
   )
   coef_seq <- matrix(0, nrow = d, ncol = length(lambda))
   results = matrix(NA, nrow = length(lambda), ncol = 5)
   # proximal gradient method
   for (s in 1:length(lambda)) {
      theta_ok <- rep(0, d) #current iteration
      theta_kk <- rep(2, d) #? previous iteration
      t_k <- 1
      k <- 0
      iter <- 0
      while (sum(abs(theta_ok - theta_kk)) > eps & (iter < max.iter)) {
         theta_k <- theta_ok
         grad <- lossgrad(
            theta = theta_k,
            beta = beta,
            A.m = A.m,
            mvc.c = mvc.c,
            y.mvc = y.mvc,
            y.is = y,
            X.is = X
         )
         proxi <- prox(
            gradL = grad,
            theta = theta_k,
            lambda = lambda[s],
            t = t_k
         )
         G <- (theta_k - proxi) / t_k
         ## line search
         if (
            loss(
               theta = theta_k - t_k * G,
               beta = beta,
               A.m = A.m,
               mvc.c = mvc.c,
               y.mvc = y.mvc,
               y.is = y,
               X.is = X
            ) >
               loss(
                  theta = theta_k,
                  beta = beta,
                  A.m = A.m,
                  mvc.c = mvc.c,
                  y.mvc = y.mvc,
                  y.is = y,
                  X.is = X
               ) -
                  t_k * t(grad) %*% G +
                  (t_k / 2) * norm(G, '2')^2
         ) {
            t_k <- t_k * tau
         } else {
            theta_ok <- proxi
            theta_kk <- theta_k
            k <- k + 1
         }
         iter <- iter+1
         if(iter==max.iter) paste0("Maximum first step iterations reached for lambda = ", lambda[s] %>% round(-log(lambda[1], 10))) %>% message()
      }
      coef_seq[, s] <- theta_ok
   }

   BIC <- apply(coef_seq, 2, loss, beta, A.m, mvc.c, y.mvc, y, X) +
      log(n) * alpha.bic * apply(coef_seq, 2, function(x) length(which(x != 0)))

   indices.on <- which(theta != 0) #only used for results table--oracle \theta does not affect method output

   return(list(
      coef = coef_seq[, which.min(BIC)],
      best.lambda = lambda[which.min(BIC)],
      # results = data.frame(lambda=lambda, ind.1 = coef_seq[indices.on[1], ], ind.2 = coef_seq[indices.on[2], ], count.on=colSums(coef_seq != 0), BIC = BIC)
      results = matrix(
         c(
            lambda %>% round(-log(lambda[1], 10)),
            coef_seq[indices.on[1], ],
            coef_seq[indices.on[2], ],
            colSums(coef_seq != 0),
            BIC
         )
      )
   ))
}
