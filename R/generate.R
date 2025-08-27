#' @export generate


generate <- function(lattice=TRUE, n, d, k=2, p.edge, s, beta=0.2, seed=0, n.burnin=30000, n.sim=100, keep.every=5,
                     verbose=FALSE){

  if(s > d) stop("Please choose s <= d.")

# source("C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/mvc.R")

  # n = 5
  # d = 3
  # p.edge = 0.30
  # s = 2
  # beta = 0.2
  # graph.seed=1
  # n.it <- 30000
  # n.sim <- 100

  require(graphon)
  require(katlabutils)
  require(igraph)

  if(lattice){
    factors <- c(1,n)
    for(p in 2:floor(sqrt(n))){
      if(!n%%p){
        factors <- c(p,n/p)
        }
      }
    if(factors[1]==1){
      n <- n+1
      message(paste0("Your choice of n is a prime number. Using n = ", n, ".") )
      for(p in 2:floor(sqrt(n))){
      if(!n%%p){
        factors <- c(p,n/p)
          }
        }
      }
    message(paste0("Using rectangular lattice with ", factors[1], " rows and ", factors[2], " columns.") )

    position.matrix <- matrix(1:n, nrow=factors[1], ncol=factors[2])

  # lattice.matrix <- matrix(NA,nrow=factors[1],ncol=factors[2])
  # for(i in 1:factors[1]){
  #   for(j in 1:factors[2]){
  #     lattice.matrix[i,j] <- paste0(i,"|",j)
  #     }
  #   }
  #
  # lattice.indices <- str_extract_all(lattice.matrix, "\\d")
  # lattice.points <- cbind(point=1:n,row=lapply(lattice.indices, '[[', 1) |> unlist() |> as.numeric(),
  #                         column=lapply(lattice.indices, '[[', 2) |> unlist() |> as.numeric())

  lattice.points <- cbind(point=1:n,row=1:factors[1],column=rep(1:factors[2],times=rep(factors[1],n/factors[1])))


  A <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    # position <- lattice.matrix[(i-1)%%factors[1] + 1,(i-1)%/%factors[1]+1] |> str_split_1('\\D') |> as.numeric()
    position <- lattice.points[i,2:3]
    distances <- abs(lattice.points[,2:3]-matrix(rep(position, n), ncol=2, byrow = TRUE) )
    neighbors <- lattice.points[rowSums(distances)==1,]
    A[i,neighbors[,1]] = 1
    }

  }



  if(!lattice){
    if(seed) set.seed(seed)
    A <- gmodel.ER(n, mode = "prob", par = p.edge, rep = 1) #reinstate once mvc algorithm is general
  }


     ## construct the covariates

  if(k < 0){
    k = -k
    message("Using uniform bound parameter k=", k, " for covariate generation.")
  }

  if(k == 0){
    set.seed(1)
    k <- rpois(1, 2)
    message("Using uniform bound parameter k=", k, " for covariate generation.")
  }

    X <- matrix(runif(n*d,min=-k,max=k),nrow=n,ncol=d)
   # X <- matrix(0, nrow=n, ncol=d)
   # sigma <- katlabutils::generate_cov_ar1(0.2, d)
   # X <- katlabutils::fast_generate_mvn(rep(0,d), sigma, n)


   # sparse vector
   theta <- rep(0, d)
   # theta[1:s] <- runif(s, 0.5, 1)*(2*rbinom(s, 1, 0.5) - 1)
   nonzero.inds <- sample(d, s)
   theta[nonzero.inds] <- sample(c(-1,1), s, replace=TRUE)
   # theta[nonzero.inds] <- sample(c(-1,1), s, replace=TRUE)*runif(s, 0.5, 1)
   predictor <- X %*% theta

      ## Gibbs Sampler
   if(keep.every%%1>0){
     keep.every <- ceiling(keep.every)
     message(paste0("keep.every must be an integer. generate() increased keep.every to ", keep.every, "."))
   }
   set.seed(Sys.time())
   # y <- rep(-1, n)
   # y <- sample(c(-1,1), n, replace=TRUE, prob=c(0.495, 0.505) )
   y <- rep(0,n)
   inds.0 <- sample(n, n/2)
   y[inds.0] <- -1
   y[-inds.0] <- 1
   m <- rep(0, n)
   for (i in 1:n) {
      m[i] <- sum(y*A[i,])
   }
   trans_pr <- rep(0, n)
   y.return <- matrix(NA, nrow=n, ncol=n.sim)

   for (t in 1:n.burnin+n.sim*keep.every) {
     # Original sampling algorithm---Gibbs sampling incorrect here
      m <- A %*% y
      trans_pr <- exp(predictor+beta*m)/(exp(predictor+beta*m)+exp(-predictor-beta*m))
      y <- 2*rbinom(n, 1, trans_pr)-1
     # indices.used <- c()
     # for(i in 1:n){
     #   vertex = sample(setdiff(1:n, indices.used), 1)
     #   indices.used <- c(indices.used,vertex)
     #   prob.vertex <- (exp(predictor+beta*m)/(exp(predictor+beta*m)+exp(-predictor-beta*m)))[vertex]
     #   y.old <- y[vertex]
     #   y[vertex] <- 2*rbinom(1,1,prob.vertex)-1
     #   # m <- A%*%y
     #   m <- m + A[,vertex]*(y[vertex]-y.old) #quicker update for m
     # }
     if(t>n.burnin & !(t-n.burnin)%%keep.every){
       iterate <- (t-n.burnin)/keep.every
       y.return[,iterate] <- y
     }

     if(!t%%1000 & verbose) print(paste0(t, " Gibbs iterations completed."))
   }

   return(list(Y=y.return, X=X, theta=theta, beta=beta, A=A, A.build=if(lattice) position.matrix else NULL, probs=trans_pr))
}
