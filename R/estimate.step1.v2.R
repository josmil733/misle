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
   theta
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
      theta_s <- rep(0, d)
      theta_s1 <- rep(100, d)
      t_s <- 1
      while (sum(abs(theta_s - theta_s1)) > eps) {
         grad <- lossgrad(
            theta = theta_s,
            beta = beta,
            A.m = A.m,
            mvc.c = mvc.c,
            y.mvc = y.mvc,
            y.is = y,
            X.is = X
         )
         proxi <- prox(
            gradL = grad,
            theta = theta_s,
            lambda = lambda[s],
            t = t_s
         )
         G <- (theta_s - proxi) / t_s
         ## line search
         if (
            loss(
               theta = theta_s - t_s * G,
               beta = beta,
               A.m = A.m,
               mvc.c = mvc.c,
               y.mvc = y.mvc,
               y.is = y,
               X.is = X
            ) >
               loss(
                  theta = theta_s,
                  beta = beta,
                  A.m = A.m,
                  mvc.c = mvc.c,
                  y.mvc = y.mvc,
                  y.is = y,
                  X.is = X
               ) -
                  t_s * t(grad) %*% G +
                  (t_s / 2) * norm(G, '2')^2
         ) {
            t_s <- t_s * tau
            theta_s1 <- theta_s
         } else {
            theta_s <- theta_s1
            theta_s1 = proxi
         }
      }
      coef_seq[, s] <- theta_s1
   }

   BIC <- apply(coef_seq, 2, loss, beta, A.m, mvc.c, y.mvc, y, X) +
      log(n) * apply(coef_seq, 2, function(x) length(which(x != 0)))

   indices.on <- which(theta != 0) #only used for results table--oracle \theta does not affect method output

   return(list(
      coef = coef_seq[, which.min(BIC)],
      best.lambda = lambda[which.min(BIC)],
      # results = data.frame(lambda=lambda, ind.1 = coef_seq[indices.on[1], ], ind.2 = coef_seq[indices.on[2], ], count.on=colSums(coef_seq != 0), BIC = BIC)
      results = matrix(
         c(
            lambda,
            coef_seq[indices.on[1], ],
            coef_seq[indices.on[2], ],
            colSums(coef_seq != 0),
            BIC
         )
      )
   ))
}
