#' @export simulate

simulate <- function(
  n,
  d,
  s,
  beta,
  k = 2,
  n.sim,
  seed = 0,
  lattice = TRUE,
  p.edge = 0.5,
  n.burnin = 5000,
  keep.every = 5,
  verbose = FALSE,
  n.lambda = 20,
  eps = .00001,
  tau = 0.8,
  sample.split = TRUE,
  p.max.iter = 6,
  compare.to.cgm = FALSE,
  optimize.cgm = FALSE,
  compare.to.vdg = FALSE,
  proposed.method = TRUE,
  inherit.data = NULL,
  r.resume = NULL,
  data.resume = NULL,
  auto.save = FALSE
) {
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
  require(listr)
  require(hdi)
  require(misle)

  # source("C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/LSW_logit2.R") #a version of Cai, Guo, and Ma's (2021) original that removes some unnecessary things and reformats code for more efficient calculation.

  # generate simulation data and parameters

  if (is.null(r.resume) & is.null(inherit.data)) {
    paste0(
      "Beginning simulation with the following parameters:
        n=",
      n,
      " | d=",
      d,
      " | s=",
      s,
      " | beta=",
      beta,
      " | n.sim=",
      n.sim
    ) |>
      message()

    data <- generate(
      lattice,
      n,
      d,
      k,
      p.edge,
      s,
      beta,
      seed,
      n.burnin,
      n.sim,
      keep.every,
      verbose
    )
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
    y.mvc <- decomp$y.ordered[1:length(mvc), ]
    y.is = decomp$y.ordered[(length(mvc) + 1):n, ]
    X.is <- decomp$X.ordered[(length(mvc) + 1):n, , ]

    # implement optional splitting of IS data into training and debiasing sets, for theoretical ease

    if (sample.split) {
      set.seed(1)
      train <- sample(
        length(mvc.c),
        ceiling(length(mvc.c) / 2),
        replace = FALSE
      )
      debias <- setdiff(1:length(mvc.c), train)
      A.m <- A[mvc, mvc.c[train]]
      dimnames(A.m) <- list(mvc, mvc.c[train])
      y.is.train <- y.is[train, ]
      X.is.train <- X.is[train, , ]
      y.is.debias <- y.is[debias, ]
      X.is.debias <- X.is[debias, , ]
    } else {
      A.m <- A.m[mvc, mvc.c] #1/29: should cause an error
    }

    to_01 <- function(x) {
      #convert y variables to {0,1} for debiasing framework
      if (all(x %in% c(-1, 1))) {
        (x + 1) / 2
      } else {
        stop("to_01() cannot accept any values except {-1,1}.")
      }
    }

    # rejection.counter <- numeric(n.sim) #0 if fail to reject H_0, 1 otherwise

    # if(proposed.method) sqe <- matrix(-1, nrow=n.sim, ncol=d)

    # if(compare.to.cgm) sqe.cgm <- matrix(-1, nrow=n.sim, ncol=d)

    # if(compare.to.vdg) sqe.vdg <- matrix(-1, nrow=n.sim, ncol=d)

    # record histogram on all nonzero parameters, plus <=2 sparse parameters
    set.seed(1)
    indices.off <- setdiff(1:d, indices.on)
    if (length(indices.off) > 0) {
      sparse.ind <- sample(
        indices.off,
        min(length(indices.off), 2),
        replace = FALSE
      )
    } else {
      sparse.ind <- integer(0)
    }

    covts.on <- (if (length(indices.on) <= 2) indices.on else indices.on[1:2])
    covts.record <- c(covts.on, sparse.ind)

    methods <- c("proposed", "CGM", "vdG")
    methods.detected <- hutils::if_else(
      c(proposed.method, compare.to.cgm, compare.to.vdg),
      methods,
      ""
    ) |>
      str_subset(".")

    # results.data = array(
    #   -1,
    #   dim = c(length(covts.record), n.sim, length(methods.detected)),
    #   dimnames = list(
    #     par.no = covts.record,
    #     sim = 1:n.sim,
    #     method = methods.detected
    #   )
    # )

    # output = list() #master output object
    # for (method in methods.detected) {}

    if (proposed.method) {
      results.hist.proposed <- tibble(
        "parameter.no" = covts.record,
        "true.val" = theta[covts.record],
        "first.step" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        "value" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        "se" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        'ci' = array(
          NA,
          dim = c(length(covts.record), n.sim, 2),
          dimnames = list(
            par.no = covts.record,
            sim = 1:n.sim,
            bound = c('lower', 'upper')
          )
        )
      )
    }

    if (compare.to.cgm) {
      results.hist.cgm <- tibble(
        "parameter.no" = covts.record,
        "true.val" = theta[covts.record],
        "first.step" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        "value" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        "se" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        'ci' = array(
          NA,
          dim = c(length(covts.record), n.sim, 2),
          dimnames = list(
            par.no = covts.record,
            sim = 1:n.sim,
            bound = c('lower', 'upper')
          )
        )
      )
    }

    if (compare.to.cgm) {
      cons = 0.05
      cons2 = 0.01 #not sure how to pick these.
      # if(d==400){cons=0.5;cons2=0.07}
      # if(d==700){cons=0.85;cons2=0.05}
      # if(d>800){cons=2;cons2=0.04}
      train.ind.cgm <- sample(1:n, ceiling(n / 2))
      debias.ind.cgm <- setdiff(1:n, train.ind.cgm)
    }

    if (compare.to.vdg) {
      results.hist.vdg <- tibble(
        "parameter.no" = covts.record,
        "true.val" = theta[covts.record],
        "value" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        "se" = matrix(NA, nrow = length(covts.record), ncol = n.sim),
        'ci' = array(
          NA,
          dim = c(length(covts.record), n.sim, 2),
          dimnames = list(
            par.no = covts.record,
            sim = 1:n.sim,
            bound = c('lower', 'upper')
          )
        )
      )
    }
  }

  if (!is.null(inherit.data)) {
    load(inherit.data)
    env_coalesce(current_env(), output)

    # sqe.cgm <- matrix(-1, nrow = n.sim, ncol = d)

    sqe.vdg <- matrix(-1, nrow = n.sim, ncol = d)

    results.hist.cgm <- tibble(
      "parameter.no" = c(indices.on, sparse.ind),
      "is.nonzero" = c(
        rep(TRUE, length(indices.on)),
        rep(FALSE, length(sparse.ind))
      ),
      "first.step" = matrix(
        NA,
        nrow = length(indices.on) + length(sparse.ind),
        ncol = n.sim
      ),
      "value" = matrix(
        NA,
        nrow = length(indices.on) + length(sparse.ind),
        ncol = n.sim
      ),
      "se" = matrix(
        NA,
        nrow = length(indices.on) + length(sparse.ind),
        ncol = n.sim
      ),
      'ci' = array(
        NA,
        dim = c(length(indices.on) + length(sparse.ind), n.sim, 2),
        dimnames = list(
          par.no = c(indices.on, sparse.ind),
          sim = 1:n.sim,
          bound = c('lower', 'upper')
        )
      )
    )

    cons = 0.05
    cons2 = 0.01 #not sure how to pick these.
    # if(d==400){cons=0.5;cons2=0.07}
    # if(d==700){cons=0.85;cons2=0.05}
    # if(d>800){cons=2;cons2=0.04}
    train.ind.cgm <- sample(1:n, ceiling(n / 2))
    debias.ind.cgm <- setdiff(1:n, train.ind.cgm)

    results.hist.vdg <- tibble(
      "parameter.no" = c(indices.on, sparse.ind),
      "is.nonzero" = c(
        rep(TRUE, length(indices.on)),
        rep(FALSE, length(sparse.ind))
      ),
      "value" = matrix(
        NA,
        nrow = length(indices.on) + length(sparse.ind),
        ncol = n.sim
      ),
      "se" = matrix(
        NA,
        nrow = length(indices.on) + length(sparse.ind),
        ncol = n.sim
      ),
      'ci' = array(
        NA,
        dim = c(length(indices.on) + length(sparse.ind), n.sim, 2),
        dimnames = list(
          par.no = c(indices.on, sparse.ind),
          sim = 1:n.sim,
          bound = c('lower', 'upper')
        )
      )
    )
  }

  interrupted <- 1 #flag to indicate whether replications stopped with error
  start <- ifelse(!is.null(r.resume), r.resume, 1)
  if (!is.null(r.resume)) {
    # attach(data.resume)
    # env_coalesce(.GlobalEnv, data.resume)
    env_coalesce(current_env(), data.resume)
    paste0(
      "Resuming simulations beginning at replication number ",
      r.resume,
      "."
    ) |>
      message()
  }

  r <- start
  while (r <= n.sim) {
    # Compare to vdg

    if (compare.to.vdg) {
      stop("compare.to.vdg temporarily disabled (10/14/2025)")
      if (r == 1) {
        fit.lasso <- hdi::lasso.proj(
          x = X[,, r],
          y = y[, 1],
          standardize = F,
          family = "binomial",
          return.Z = TRUE
        )
        Z <- fit.lasso$Z
      }
      fit.lasso <- hdi::lasso.proj(
        x = X,
        y = y[, r],
        Z = Z,
        standardize = F,
        family = "binomial",
        return.Z = FALSE
      )
      theta.tilde.vdg <- fit.lasso$bhat
      sqe.vdg[r, ] <- (theta.tilde.vdg - theta)^2
      se.vdg <- fit.lasso$se
      for (j in covts.record) {
        if (j %in% results.hist.vdg$parameter.no) {
          results.hist.vdg$value[
            which(results.hist.vdg$parameter.no == j),
            r
          ] = theta.tilde.vdg[j]
          results.hist.vdg$se[
            which(results.hist.vdg$parameter.no == j),
            r
          ] = se.vdg[j]
          results.hist.vdg$ci[
            which(results.hist.vdg$parameter.no == j),
            r,
          ] = c(
            theta.tilde.vdg[j] - 1.96 * se.vdg[j],
            theta.tilde.vdg[j] + 1.96 * se.vdg[j]
          )
        }
      }
    }

    # Step 1: construct \ell_1-penalized MLE

    if (proposed.method) {
      theta.hat <- estimate.step1(
        grid.lambda = list(from = 0.01, to = 0.1, length.out = 20),
        d = d,
        eps = eps,
        y.mvc = y.mvc[, r],
        y = if (sample.split) y.is.train[, r] else y.is[, r],
        X = if (sample.split) X.is.train[,, r] else X.is[,, r],
        beta = beta,
        A.m = A.m,
        mvc.c = mvc.c,
        tau = 0.8
      )
    }

    if (compare.to.cgm) {
      y.cgm.train <- to_01(y[train.ind.cgm, r])
      y.cgm.debias <- to_01(y[debias.ind.cgm, r])
      X.cgm.train <- 2 * X[train.ind.cgm, , r]
      X.cgm.debias <- 2 * X[debias.ind.cgm, , r]
      if (optimize.cgm) {
        theta.hat.cgm <- cgm.inference1(
          X = X.cgm.train,
          y = y.cgm.train,
          lambda = NULL
        )
      } else {
        # theta.hat.cgm <- cgm.inference1(X = X.cgm.train, y = to_01(y.cgm.train), lambda = cons2*sqrt(log(d)/n))
        # theta.hat.cgm <- cgm.inference1(X = X.cgm.train, y = y.cgm.train, lambda = cons2*sqrt(2*log(d)/n))
        theta.hat.cgm <- glmnet(
          X.cgm.train,
          y.cgm.train,
          family = "binomial",
          alpha = 1,
          intercept = FALSE,
          lambda = cons2 * sqrt(2 * log(d) / n),
          standardize = FALSE
        )$beta |>
          as.vector()
      }
    }

    # Step 2: debias theta.hat

    if (proposed.method) {
      theta.tilde <- rep(NA, d)
      se <- rep(NA, d)
    }

    if (compare.to.cgm) {
      theta.tilde.cgm <- rep(NA, d)
      se.cgm <- rep(NA, d)
    }

    error.flag <- 0

    for (j in covts.record) {
      if (proposed.method) {
        # debias <- estimate.step2(theta.hat=theta.hat, beta = beta, A.m = A.m,
        # X=if(sample.split) X.is.train else X.is,
        # y=if(sample.split) to_01(y.is.train[,r]) else to_01(y.is[,r]),
        # y.mvc = y.mvc[,r], predictor=j)

        debias <- estimate.step2(
          theta.hat = theta.hat,
          beta = beta,
          A.m = A.m,
          X = if (sample.split) X.is.debias[,, r] else X.is[,, r],
          y = if (sample.split) to_01(y.is.debias[, r]) else to_01(y.is[, r]),
          y.mvc = y.mvc[, r],
          predictor = j,
          p.max.iter = p.max.iter
        )

        if ('error' %in% names(debias)) {
          error.flag <- 1
          err.source <- 'MISLE'
        }
      }

      if (compare.to.cgm) {
        debias.cgm <- cgm.inference2(
          theta.hat.cgm,
          X.cgm.debias,
          y.cgm.debias,
          j
        )
        #   if('error' %in% names(debias.cgm)){
        #     error.flag <- 1
        #     if('error' %in% names(debias) & !'error' %in% names(debias.cgm)) err.source='MISLE'
        #     if('error' %in% names(debias.cgm) & !'error' %in% names(debias)) err.source='CGM'
        #     if('error' %in% names(debias) & 'error' %in% names(debias.cgm)) err.source='MISLE and CGM'
        # }
        if ('error' %in% names(debias.cgm)) {
          error.flag <- 1
          err.source <- 'CGM'
        }
      }

      if (error.flag) {
        wrn <- paste0(
          'Error detected in direction vector computation for ',
          err.source,
          ' for replicate r=',
          r,
          ' and covariate j=',
          j,
          '. Generating new data and attempting this iteration again.'
        )
        message(wrn)
        break
      }

      if (proposed.method) {
        theta.tilde[j] <- debias$theta.tilde
        se[j] <- debias$se
        if (j %in% results.hist.proposed$parameter.no) {
          results.hist.proposed$first.step[
            which(results.hist.proposed$parameter.no == j),
            r
          ] = theta.hat[j]
          results.hist.proposed$value[
            which(results.hist.proposed$parameter.no == j),
            r
          ] = theta.tilde[j]
          results.hist.proposed$se[
            which(results.hist.proposed$parameter.no == j),
            r
          ] = se[j]
          results.hist.proposed$ci[
            which(results.hist.proposed$parameter.no == j),
            r,
          ] = c(
            theta.tilde[j] - 1.96 * se[j],
            theta.tilde[j] + 1.96 * se[j]
          )
        }
      }

      if (compare.to.cgm) {
        theta.tilde.cgm[j] <- debias.cgm$theta.tilde
        se.cgm[j] <- debias.cgm$se
        if (j %in% results.hist.cgm$parameter.no) {
          results.hist.cgm$first.step[
            which(results.hist.cgm$parameter.no == j),
            r
          ] = theta.hat.cgm[j]
          results.hist.cgm$value[
            which(results.hist.cgm$parameter.no == j),
            r
          ] = theta.tilde.cgm[j]
          results.hist.cgm$se[
            which(results.hist.cgm$parameter.no == j),
            r
          ] = se.cgm[j]
          results.hist.cgm$ci[
            which(results.hist.cgm$parameter.no == j),
            r,
          ] = c(
            theta.tilde.cgm[j] - 1.96 * se.cgm[j],
            theta.tilde.cgm[j] + 1.96 * se.cgm[j]
          )
        }
      }

      # if(proposed.method){
      #   p.values <- 2*pnorm(-abs(sqrt(n)*theta.tilde/se))
      #   # debias <- deb.lasso(x=X.is, y=y.is, lasso_est=theta.hat, inference=TRUE)
      #   # rejection.counter[r] <- 1*(any(debias$p.value <= 0.05/d)) #Bonferroni correction
      #   rejection.counter[r] <- 1*(any(p.values <= 0.05/d)) #Bonferroni correction
      #   sqe[r,] <- (theta.tilde-theta)^2
      # }
      #
      # if(compare.to.cgm) sqe.cgm[r,] <- (theta.tilde.cgm-theta)^2

      # if(verbose & !j%%max(1,floor(d/5)) & (proposed.method | compare.to.cgm)) message( paste0("replication ", r, ": ", j, " covariates complete."))
      # }
    }

    if (error.flag) {
      data.regen <- regenerate(r, y[, n.sim])
      data$Y[, r] = data.regen$y
      data$X[,, r] = data.regen$X
      y.ordered <- numeric(n)
      y.ordered[1:length(mvc)] <- data.regen$y[mvc]
      y.ordered[(length(mvc) + 1):n] <- data.regen$y[mvc.c]
      names(y.ordered) <- c(mvc, mvc.c)
      y.mvc[, r] <- y.ordered[1:length(mvc)]
      y.is[, r] = y.ordered[(length(mvc) + 1):n]
      X.ordered <- matrix(NA, nrow = n, ncol = d)
      X.ordered[1:length(mvc), ] <- data.regen$X[mvc, ]
      X.ordered[(length(mvc) + 1):n, ] <- data.regen$X[mvc.c, ]
      rownames(X.ordered) <- c(mvc, mvc.c)
      X.is[,, r] = X.ordered[(length(mvc) + 1):n, ]
      if (sample.split) {
        y.is.train[, r] <- y.is[train, r]
        y.is.debias[, r] <- y.is[setdiff(1:length(mvc.c), train), r]
        X.is.train[,, r] <- X.is[train, , r]
        X.is.debias[,, r] <- X.is[setdiff(1:length(mvc.c), train), , r]
      }
      next
    }

    if (auto.save & (proposed.method | compare.to.cgm)) {
      if (!r %% min(5, ceiling(n.sim / 5))) {
        results <- current_env()
        file.name <- paste0(
          "n",
          n,
          "-d",
          d,
          "-beta",
          beta,
          "-s",
          s,
          "---partial-r=",
          r,
          ".RData"
        )
        # save(results, file=paste0('C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/Results/', file.name) )
        save(results, file = paste0(getwd(), "/", file.name))
        if (r == min(5, ceiling(n.sim / 5))) {
          paste0(
            "Saving simulation data to ",
            getwd(),
            ". Temporary data files will be named in the format n",
            n,
            "-d",
            d,
            "-beta",
            beta,
            "-s",
            s,
            "---partial-r=.RData"
          ) |>
            message()
        } else {
          message("Save successful after ", r, " replications.")
        }
      }
    }

    if (verbose & !r %% 5) {
      message(paste0("replication ", r, " complete."))
    }

    r = r + 1
    start <- r #for ease of debugging

    ## Use (e.g.) Javanmard and Montanari 2014 as starting point for our debiasing approach

    # asymp.var <- asymptotic variance of theta.tilde
    # test.stat <- sqrt(n)*theta.tilde/sqrt(asymp.var) #under H_0: \theta=0 (2-sided)
    # p.value <- 2*pnorm(test.stat)
    # rejection.counter[r] <- 1*(p.value <= 0.05)
  }

  interrupted = 0
  return(current_env())
}
