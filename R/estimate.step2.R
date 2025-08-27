#' @export f
#' @export fp
#' @export w
#' @export getmode
#' @export estimate.step2

f = function(x){
   exp(x)/(exp(x) + exp(-x))
}
fp = function(x){
   2*f(x)/(1+exp(2*x))
}
w = function(x){
   fp(x)/( f(x)*(1-f(x)) )
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



estimate.step2 <- function(theta.hat, beta, A.m, X, y.mvc, y, predictor, p.max.iter){
    ### Option 1: search tuning parameter with steps determined by the ill conditioned case (n=p/2)
  ### Option 2: search tuning parameter with maximum 10 steps.
  ####### Option 3: fixed tuning parameter and this is not recommended without exploring the tuning parameter selection

  d <- ncol(X)

  if(!(predictor %in% (1:d))) stop("predictor should be an integer between 1 and d.")

  # y <- y.is.train[,r]
  # X <- X.is.train
  mu=NULL
  step=NULL
  resol=1.5
  cons=2.01
  # maxiter=6

  n <- length(y) #2/13: should this be changed?

    loading = rep(0,d)
    loading[predictor]=1
    lin.pred.hat= as.vector(2*beta*t(A.m)%*%y.mvc) + X%*%(theta.hat)
    X.weight = diag(c(sqrt(fp(lin.pred.hat)*w(lin.pred.hat)))) %*% X
  #####################################################################################################
  ################## Correction step
   if ((n>=6*d)){
      sigma.hat <- (1/n)*(t(X.weight)%*%X.weight);
      tmp <- eigen(sigma.hat)
      tmp <- min(tmp$values)/max(tmp$values)
   }else{
      tmp <- 0
   }

   if ((n>=6*d)&&(tmp>=1e-4)){
      direction <- solve(sigma.hat)%*%loading
   }else{
      if(n>0.5*d){
        # print(paste0('n>0.5d called.'))
         ### for option 1
         if(is.null(step)){
            step.vec<-rep(NA,3)
            for(t in 1:3){
               index.sel<-sample(1:n,size=ceiling(0.5*min(n,d)), replace=FALSE) #sample smaller sample

               Direction.Est.temp<-Direction_searchtuning(Xc=X.weight[index.sel,],loading=loading,mu=NULL, resol=resol, maxiter=p.max.iter)
               if(Direction.Est.temp$error.code){
                  return(list(error=1))
               }
               step.vec[t]<-Direction.Est.temp$step
            }
            step<-getmode(step.vec)
         }
         # print(paste("step is", step))
         Direction.Est<-Direction_fixedtuning(X.weight,loading,mu=sqrt(cons*log(d)/n)*resol^{-(step-1)})
        if(Direction.Est$error.code){
          return(list(error=1))
          }
      }else{
         ### for option 2
        # print(paste0("n>0.5d not called."))
         Direction.Est<-Direction_searchtuning(X.weight,loading,mu=NULL, resol, maxiter=p.max.iter)
          if(Direction.Est$error.code){
            return(list(error=1))
            }
         step<-Direction.Est$step
         # print(paste("step is", step))
      }
      direction<-Direction.Est$proj
   }

   weighed.residual=(y - f(lin.pred.hat))*w(lin.pred.hat)

   correction = sum((X%*%direction)*weighed.residual)/n;
   debias.est=theta.hat[predictor]+correction

   X.weight2 = diag(c(w(lin.pred.hat)*sqrt(f(lin.pred.hat)*(1-f(lin.pred.hat))))) %*% X
   # se<-sqrt(mean((X.weight2%*%direction)^2))/sqrt(n)
   se<-(mean( (X.weight2%*%direction)^2 )/n) |> sqrt() #changed 8/4/2025
   #sd
   #mean((y - exp(lin.pred.hat)/(1+ exp(lin.pred.hat)))^2)
   #mean((y - exp(exp_val)/(1+ exp(exp_val)))^2)
   #mean(exp(lin.pred.hat)/(1+ exp(lin.pred.hat))^2)
   #c(linear.plugin+correct-1.96*sd,linear.plugin+correct+1.96*sd)
   returnList <- list("theta.tilde" = debias.est,
                      "j"=predictor,
                      "se" = se,
                      "u.hat"=direction,
                      "step"=step,
                      "theta.hat"=theta.hat
   )
   return(returnList)
}
