## loss function in ISL Ising Model
loss <- function(theta, beta, A.m, mvc.c, y.mvc, y.is, X.is){
   v <- (2*beta*t(A.m)%*%y.mvc + X.is%*%theta) |> as.vector()
   names(v) <- colnames(X.is)
   -(1/length(mvc.c))*sum(y.is*v - log(cosh(v)))
}

## loss gradient of ISL Ising Model
lossgrad <- function(theta, beta, A.m, mvc.c, y.mvc, y.is, X.is){
   v <- (2*beta*t(A.m)%*%y.mvc + X.is%*%theta) |> as.vector()
   names(v) <- colnames(X.is)
   (1/length(mvc.c))*colSums(tanh(v)*X.is-y.is*X.is)
}

## proximal gradient descent
prox <- function(gradL, theta, lambda = lambda, t){
      b <- rep(0, length(theta))
      for (j in 1:length(theta)) {
         if ((theta[j]-t*gradL[j])>t*lambda){
            b[j] <- theta[j]-t*gradL[j] - lambda*t
         }
         if ((theta[j]-t*gradL[j])<t*lambda){
            b[j] <- theta[j]-t*gradL[j] + lambda*t
         }
         if (abs(theta[j]-t*gradL[j])<=lambda*t){
            b[j] <- 0
         }
      }
      return(b)
   }
