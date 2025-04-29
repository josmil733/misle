#' @export Direction_fixedtuning
#' @export Direction_searchtuning
#' @export cgm.inference1
#' @export cgm.inference2


library(MASS)
library(CVXR)
library(AER)
library(Matrix)
library(glmnet)



f = function(x){
  exp(x)/(1+exp(x))
}
fp = function(x){
  exp(x)/(1+exp(x))^2
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Direction_fixedtuning<-function(Xc,loading, mu=NULL){
  pp<-ncol(Xc)
  n<-nrow(Xc)
  if(is.null(mu)){
    mu<-sqrt(2.01*log(pp)/n)
  }
  loading.norm<-sqrt(sum(loading^2))
  H<-cbind(loading/loading.norm,diag(1,pp))
  v<-Variable(pp+1)
  obj<-1/4*sum((Xc%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob<-Problem(Minimize(obj))
<<<<<<< HEAD
<<<<<<< HEAD
  result<-CVXR::solve(prob)
=======
  result<-solve(prob)
>>>>>>> db613dc (restore original direction solvers)
=======
  result<-solve(prob)
>>>>>>> db613dc (restore original direction solvers)
  print("fixed mu")
  print(mu)
  #print(result$value)
  opt.sol<-result$getValue(v)
  cvxr_status<-result$status
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  returnList <- list("proj"=direction)
  return(returnList)
}



Direction_searchtuning<-function(Xc,loading,mu=NULL, resol, maxiter){
  pp<-ncol(Xc)
  n<-nrow(Xc)
  tryno = 1;
  opt.sol = rep(0,pp);
  lamstop = 0;
  cvxr_status = "optimal";

  mu = sqrt(2.01*log(pp)/n);
  #mu.initial= mu;
  while (lamstop == 0 && tryno < maxiter){
    ###### This iteration is to find a good tuning parameter
    #print(mu);
    lastv = opt.sol;
    lastresp = cvxr_status;
    loading.norm<-sqrt(sum(loading^2))
    H<-cbind(loading/loading.norm,diag(1,pp))
    v<-Variable(pp+1)
    obj<-1/4*sum((Xc%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    prob<-Problem(Minimize(obj))
    result<-solve(prob)
    #print(result$value)
    opt.sol<-result$getValue(v)
    cvxr_status<-result$status
    #print(cvxr_status)
    if(tryno==1){
      if(cvxr_status=="optimal"){
        incr = 0;
        mu=mu/resol;
        temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        initial.sd<-sqrt(sum((Xc%*% temp.vec)^2)/(n)^2)*loading.norm
        temp.sd<-initial.sd
      }else{
        incr = 1;
        mu=mu*resol;
      }
    }else{
      if(incr == 1){ ### if the tuning parameter is increased in the last step
        if(cvxr_status=="optimal"){
          lamstop = 1;
        }else{
          mu=mu*resol;
        }
      }else{
        if(cvxr_status=="optimal"&&temp.sd<3*initial.sd){
          mu = mu/resol;
          temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          temp.sd<-sqrt(sum((Xc%*% temp.vec)^2)/(n)^2)*loading.norm
          #print(temp.sd)
        }else{
          mu=mu*resol;
          opt.sol=lastv;
          lamstop=1;
          tryno=tryno-1
        }
      }
    }
    tryno = tryno + 1;
  }
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step<-tryno-1
  print(step)
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}



cgm.inference1<-function(X,y,lambda=NULL){

  X<-as.matrix(X)

  ### implement logistic Lasso

  fit = glmnet(X, y,  family = "binomial", alpha = 1, intercept=FALSE,
               lambda = lambda, standardize=F)
  return(theta.hat=as.vector(fit$beta))
}

cgm.inference2<-function(theta.hat, X,y,predictor){
  ### Option 1: search tuning parameter with steps determined by the ill conditioned case (n=p/2)
  ### Option 2: search tuning parameter with maximum 10 steps.
  ####### Option 3: fixed tuning parameter and this is not recommended without exploring the tuning parameter selection

  d <- ncol(X)

  if(!(predictor %in% (1:d))) stop("predictor should be an integer between 1 and d.")

  mu=NULL
  step=NULL
  resol=1.5
  cons=2.01
  maxiter=6

  n <- length(y) #2/13: should this be changed?

    loading = rep(0,d)
    loading[predictor]=1
    lin.pred.hat=X%*%(theta.hat)
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
      ### for option 1
      if(is.null(step)){
        step.vec<-rep(NA,3)
        for(t in 1:3){
          index.sel<-sample(1:n,size=ceiling(0.5*min(n,d)), replace=FALSE) #sample smaller sample

          Direction.Est.temp<-Direction_searchtuning(X.weight[index.sel,],loading,mu=NULL, resol, maxiter)
          step.vec[t]<-Direction.Est.temp$step
        }
        step<-getmode(step.vec)
      }
      # print(paste("step is", step))
      Direction.Est<-Direction_fixedtuning(X.weight,loading,mu=sqrt(cons*log(d)/n)*resol^{-(step-1)})
    }else{
      ### for option 2
      Direction.Est<-Direction_searchtuning(X.weight,loading,mu=NULL, resol, maxiter)
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
   se<-mean( (X.weight2%*%direction)^2 ) |> sqrt()



  returnList <- list("theta.tilde" = debias.est,
                     "j" = predictor,
                     "se" = se,
                     "u.hat"=direction,
                     "step"=step,
                     "theta.hat"=theta.hat
  )
  return(returnList)
}
