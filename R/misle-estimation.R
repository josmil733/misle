#' @export simulate

simulate <- function(n, d, s, beta, n.sim, seed=0, lattice=TRUE,  p.edge = 0.5, n.burnin=30000, keep.every=5, verbose=TRUE,
            n.lambda=20, eps = .00001, tau=0.8, sample.split=TRUE,
            compare.to.cgm=TRUE,
            r.resume=NULL, data.resume=NULL){


# simulate(lattice=TRUE, n=20, d=2, p.edge=0.5, s=1, beta=0.05, seed=0, n.burnin=2000, n.sim=2, keep.every=5, verbose=TRUE, n.lambda=20, eps=1e-5, tau=0.8, sample.split=TRUE)

require(MASS)
require(glmnet)
require(devtools)
require(katlabutils)
require(graphon)
require(CVXR)
require(tidyverse)
require(tibble)
require(rmarkdown)
require(rlang)
require(misle)

# source("C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/LSW_logit2.R") #a version of Cai, Guo, and Ma's (2021) original that removes some unnecessary things and reformats code for more efficient calculation.

    # generate simulation data and parameters

if(is.null(r.resume)){

  paste0("Beginning simulation with the following parameters:
        n=",n," | d=",d," | s=",s," | beta=",beta," | n.sim=",n.sim) |> message()

  data <- generate(lattice,n,d,p.edge,s,beta,seed,n.burnin,n.sim,keep.every,verbose)
  n <- nrow(data$X)
  d <- ncol(data$X)
  A <- data$A
  A.build <- data$A.build
  X <- data$X
  y <- data$Y
  theta <- data$theta
  beta <- data$beta
  indices.on <- which(theta != 0)



  # decompose A,y based on minimum vertex cover computation

  decomp <- misle::decompose(A, A.build, y, X, lattice)
  mvc <- decomp$mvc
  mvc.c <- setdiff(1:n, mvc)
  # A.m <- A[mvc, mvc.c]
  y.mvc <- decomp$y.ordered[1:length(mvc),]
  y.is = decomp$y.ordered[(length(mvc)+1):n,]
  X.is <- decomp$X.ordered[(length(mvc)+1):n,]



  # implement optional splitting of IS data into training and debiasing sets, for theoretical ease

  if(sample.split){
       set.seed(1)
       train <- sample(length(mvc.c), ceiling(length(mvc.c)/2), replace=FALSE)
       debias <- setdiff(1:length(mvc.c), train)
       A.m <- A[mvc, mvc.c[train]]
       dimnames(A.m) <- list(mvc, mvc.c[train])
       y.is.train <- y.is[train,]
       X.is.train <- X.is[train,]
       y.is.debias <- y.is[debias,]
       X.is.debias <- X.is[debias,]
  } else {
    A.m <- A.m[mvc, mvc.c] #1/29: should cause an error
  }

      to_01 <- function(x){ #convert y variables to {0,1} for debiasing framework
        if(all( x %in% c(-1,1) ) ) (x+1)/2 else stop("to_01() cannot accept any values except {-1,1}.")
      }

      rejection.counter <- numeric(n.sim) #0 if fail to reject H_0, 1 otherwise

      sqe <- matrix(-1, nrow=n.sim, ncol=d)

      if(compare.to.cgm) sqe.cgm <- matrix(-1, nrow=n.sim, ncol=d)

      # record histogram on all nonzero parameters, plus <=2 sparse parameters
      set.seed(1)
      indices.off <- setdiff(1:d, indices.on)
      if(length(indices.off)>0){
        sparse.ind <- sample(indices.off, min(length(indices.off), 2), replace=FALSE)
      } else sparse.ind <- integer(0)

      results.hist <- tibble("parameter.no" = c(indices.on,sparse.ind), is.nonzero = c(rep(TRUE, length(indices.on)), rep(FALSE, length(sparse.ind))), value = matrix(NA, nrow=length(indices.on)+length(sparse.ind), ncol=n.sim),
                                 se = matrix(NA, nrow=length(indices.on)+length(sparse.ind), ncol=n.sim))

      if(compare.to.cgm) results.hist.cgm <- tibble("parameter.no" = c(indices.on,sparse.ind), is.nonzero = c(rep(TRUE, length(indices.on)), rep(FALSE, length(sparse.ind))), value = matrix(NA, nrow=length(indices.on)+length(sparse.ind), ncol=n.sim),
                                 se = matrix(NA, nrow=length(indices.on)+length(sparse.ind), ncol=n.sim))

}

    interrupted <- 1 #flag to indicate whether replications stopped with error
    start <- ifelse(!is.null(r.resume), r.resume, 1)
    if(!is.null(r.resume)){
      # attach(data.resume)
      env_coalesce(.GlobalEnv, data.resume)
      paste0("Resuming simulations beginning at replication number ", r.resume, ".") |> message()
    }

   for (r in start:n.sim) {

      # Step 1: construct \ell_1-penalized MLE

      theta.hat <- estimate.step1(grid.lambda=list(from=0.01, to=0.1, length.out=20),
                                  d=d, eps=eps, y.mvc = y.mvc[,r], y=if(sample.split) y.is.train[,r] else y.is[,r],
                           X=if(sample.split) X.is.train else X.is,
                           beta = beta, A.m = A.m, mvc.c = mvc.c,
                           tau=0.8)

      # Step 2: debias theta.hat

      theta.tilde <- numeric(d)
      se <- numeric(d)

      for(j in 1:d){
        debias <- estimate.step2(theta.hat=theta.hat, beta = beta, A.m = A.m,
                      X=if(sample.split) X.is.train else X.is,
                      y=if(sample.split) to_01(y.is.train[,r]) else to_01(y.is[,r]),
                      y.mvc = y.mvc[,r], predictor=j)

        theta.tilde[j] <- debias$theta.tilde
        se[j] <- debias$se #12/4: may need to update this calculation.
        if(j %in% results.hist$parameter.no){
          results.hist$value[which(results.hist$parameter.no==j), r] = theta.tilde[j]
          results.hist$se[which(results.hist$parameter.no==j), r] = se[j]
        }
        if(verbose & !j%%max(1,floor(d/5)) ) print( paste0("replication ", r, ": ", j, " covariates complete."))
      }

      p.values <- 2*pnorm(-abs(sqrt(n)*theta.tilde/se))

      # debias <- deb.lasso(x=X.is, y=y.is, lasso_est=theta.hat, inference=TRUE)
      # rejection.counter[r] <- 1*(any(debias$p.value <= 0.05/d)) #Bonferroni correction
      rejection.counter[r] <- 1*(any(p.values <= 0.05/d)) #Bonferroni correction

      sqe[r,] <- (theta.tilde-theta)^2



      if(compare.to.cgm){

        cons=0.05; cons2=0.01 #not sure how to pick these.
        # if(d==400){cons=0.5;cons2=0.07}
        # if(d==700){cons=0.85;cons2=0.05}
        # if(d>800){cons=2;cons2=0.04}

        train.cgm <- sample(1:n, ceiling(n/2))
        debias.cgm <- setdiff(1:n, train.cgm)
        y.cgm.train <- y[train.cgm, r]
        y.cgm.debias <- y[debias.cgm, r]
        X.cgm.train <- X[train.cgm,]
        X.cgm.debias <- X[debias.cgm,]

        theta.tilde.cgm <- numeric(d)
        se.cgm <- numeric(d)

        theta.hat.cgm <- cgm.inference1(X = X.cgm.train, y = to_01(y.cgm.train), lambda = cons2*sqrt(log(d)/n))

        for(j in 1:d){
          # cgm.fit <- cgm.mod.jm(X=X, y=y[,r], nlambda=n.lambda, lambda.min.ratio=0.01, intercept=FALSE,
                                      # mu=NULL, step=NULL, resol=1.5, cons=2.01, maxiter=6)
          debias.cgm <- cgm.inference2(theta.hat, X.cgm.debias, y.cgm.debias, j)
          theta.tilde.cgm[j] <- debias.cgm$theta.tilde
          se.cgm[j] <- debias.cgm$se

        if(j %in% results.hist.cgm$parameter.no){
          results.hist.cgm$value[which(results.hist.cgm$parameter.no==j), r] = theta.tilde.cgm[j]
          results.hist.cgm$se[which(results.hist.cgm$parameter.no==j), r] = se.cgm[j]
        }
          if(verbose & !j%%max(1,floor(d/5)) ) print( paste0("replication ", r, ": ", j, " covariates complete (CGM)."))
        }

        sqe.cgm[r,] <- (theta.tilde.cgm-data$theta)^2

      }


  if(!r%%min(5, ceiling(n.sim/5))){
    results <- current_env()
    file.name <- paste0("n",n,"-d",d,"-beta",beta,"-s",s,"---partial-r=",r,".RData")
    # save(results, file=paste0('C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/Results/', file.name) )
    save(results, file=paste0(getwd(), "/", file.name) )
    if(r==min(5, ceiling(n.sim/5))){
      paste0("Saving simulation data to ", getwd(), ". Temporary data files will be named in the format n",n,"-d",d,"-beta",beta,"-s",s,"---partial-r=.RData") |> message()
    } else {
    paste0("Save successful after ", r, " replications.") |> message()
      }

  }


      ## Use (e.g.) Javanmard and Montanari 2014 as starting point for our debiasing approach

      # asymp.var <- asymptotic variance of theta.tilde
      # test.stat <- sqrt(n)*theta.tilde/sqrt(asymp.var) #under H_0: \theta=0 (2-sided)
      # p.value <- 2*pnorm(test.stat)
      # rejection.counter[r] <- 1*(p.value <= 0.05)
   }



  interrupted = 0
  return(current_env())
}
