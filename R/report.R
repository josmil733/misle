#' @export report

report <- function(n, d, s, beta, n.sim, p.edge=NULL, seed=0, lattice=TRUE, n.burnin=30000, keep.every=5, verbose=TRUE,
            n.lambda=20, eps = .00001, tau=0.8, sample.split=TRUE,
            compare.to.cgm=FALSE, optimize.cgm=TRUE, compare.to.vdg=TRUE, proposed.method=TRUE, inherit.data=NULL,
            simulation.path=NULL, note=NULL, ec=FALSE, results.dir=getwd()
            ){

 # perform the simulations with varying n (=400, ... 1000), covariate dimensions (take d= n/10), sparsity level (s=2,5,10). Change beta around some as well. Repeat experiments B=50 times and report what proportion of confidence intervals cover the true parameters. Write a description along with the plots of all the results you received

require(rlang)
require(stringr)
require(rmarkdown)
require(tidyverse)
require(hutils)

if(!is.null(simulation.path)){
  args.keep <- c('simulation.path', 'note' ,'ec', 'results.dir')
  rm(list=setdiff(ls(current_env()), args.keep), pos=current_env())
  sim.active <- FALSE
  load(simulation.path)
  if('results' %in% ls()){ #'results' is the name of an environment containing an incomplete simulation
  # if('results' %in% ls(output)){ #'results' is the name of an environment containing an incomplete simulation

    # env_names(results) |> as.list() %>% lapply(env_get, results)

    r.resume <- env_get(results, 'r')+1
    if(r.resume <= env_get(results, 'n.sim') ){
      env_coalesce(current_env(), results)
      # env_coalesce(.GlobalEnv, results)
      # env_coalesce(caller_env(3))
      output <- simulate(n=n, d=d, s=s, beta=beta, n.sim=n.sim, p.edge=p.edge, seed=seed, n.burnin=n.burnin, keep.every=keep.every, verbose=verbose,
           n.lambda=n.lambda, eps=eps, tau=tau, sample.split=sample.split, compare.to.cgm=compare.to.cgm, optimize.cgm=optimize.cgm, compare.to.vdg=compare.to.vdg, proposed.method=proposed.method, inherit.data=inherit.data,
           r.resume=r.resume, data.resume=results)
      sim.active = TRUE
        } else {output <- results}
      # } else {output <- current_env()}
  } else {env_bind(output, simulation.path=simulation.path)}
    } else {
      if(!is.null(inherit.data)){
        if(!str_detect(inherit.data, '.RData$')){
          stop("inherit.data argument must be a .RData file containing data to use in the simulation.")
        }
        load(inherit.data)
        env_coalesce(current_env(), output)
        message('Inheriting the data used for the simulations in ', inherit.data, ".")
        rm(output)
      }
output <- simulate(n=n, d=d, s=s, beta=beta, n.sim=n.sim, p.edge=p.edge, seed=seed, n.burnin=n.burnin, keep.every=keep.every, verbose=verbose,
           n.lambda=n.lambda, eps=eps, tau=tau, sample.split=sample.split, compare.to.cgm=compare.to.cgm, optimize.cgm=optimize.cgm, compare.to.vdg=compare.to.vdg, proposed.method=proposed.method, inherit.data=inherit.data)
  sim.active=TRUE
  }

  env_coalesce(current_env(), output)

  directory <- paste0(results.dir, "/")
  methods <- c("proposed", "CGM", "vdG")
  methods.detected <- hutils::if_else(c(proposed.method, compare.to.cgm, compare.to.vdg),
                                      methods,
                                      "") |> str_subset(".")

  # file.param.name <- paste0('n',n, '-d',d, '-beta',beta, '-s', s)
  file.param.name <- paste0('n',env_get(output, 'n'), '-d',env_get(output, 'd'), '-beta',env_get(output, 'beta'), '-s',env_get(output, 's') )
  file.method.name <- paste0(methods.detected, collapse="+")
  file.name <- paste0(file.param.name, "---", file.method.name)
  folder <- paste0(directory, ifelse(ec, 'EC---', ''), file.name)
  setwd(directory)

  if(file.exists(folder) & sim.active){
    similar.names <- list.dirs(recursive=FALSE) |> str_subset(file.name) |> str_remove("\\./")
    if(length(similar.names)>1){
      latest.version <- similar.names |> str_subset("v\\d") |> str_extract("\\d$") |> as.numeric() |> max()
    } else {latest.version <- 1}
    location <- paste0(file.name, '---v', latest.version+1)
    folder=paste0(folder, '---v', latest.version+1)
    paste0("Simulation results for these parameters and methods already exist. Saving results to folder ", location, ".") |> message()
  }

  if(sim.active){
    dir.create(folder)
  }
  env_bind(output, date=Sys.Date(),
              simulation.path=simulation.path, note=note, ec=ec, folder=folder, file.param.name=file.param.name)
  save(output, file=paste0(folder, '/sim-data.RData'))
  render('../results-template.Rmd', output_file=paste0(folder, '/report.pdf'), envir=output)
}
