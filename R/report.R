#' @export report

report <- function(n, d, s, beta, n.sim, p.edge=NULL, seed=0, lattice=TRUE, n.burnin=30000, keep.every=5, verbose=TRUE,
            n.lambda=20, eps = .00001, tau=0.8, sample.split=TRUE,
            compare.to.cgm=TRUE,
            simulation.path=NULL, note=NULL, ec=FALSE
            ){

 # perform the simulations with varying n (=400, ... 1000), covariate dimensions (take d= n/10), sparsity level (s=2,5,10). Change beta around some as well. Repeat experiments B=50 times and report what proportion of confidence intervals cover the true parameters. Write a description along with the plots of all the results you received

require(rlang)
require(stringr)

if(!is.null(simulation.path)){
  load(simulation.path)
  if('results' %in% ls()){ #'results' is the name of an environment containing an incomplete simulation
    r.resume <- env_get(results, 'r')+1
    if(r.resume <= env_get(results, 'n.sim') ){
      output <- simulate(n=n, d=d, s=s, beta=beta, n.sim=n.sim, p.edge=p.edge, seed=seed, n.burnin=n.burnin, keep.every=keep.every, verbose=verbose,
           n.lambda=n.lambda, eps=eps, tau=tau, sample.split=sample.split, compare.to.cgm=compare.to.cgm, r.resume=r.resume, data.resume=results)
        } else {output <- results}
      } else {output <- current_env()}
    } else {
output <- simulate(n=n, d=d, s=s, beta=beta, n.sim=n.sim, p.edge=p.edge, seed=seed, n.burnin=n.burnin, keep.every=keep.every, verbose=verbose,
           n.lambda=n.lambda, eps=eps, tau=tau, sample.split=sample.split, compare.to.cgm=compare.to.cgm)
  }



  directory <- "C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/Results/"
  file.name.base <- paste0('n',env_get(output, 'n'), '-d',env_get(output,'d'), '-beta',env_get(output,"beta"), '-s', env_get(output,'s'))
  folder <- paste0(directory, ifelse(ec, 'EC---', ''), file.name.base)
  setwd(directory)

  if(file.exists(folder)){
    similar.names <- list.dirs(recursive=FALSE) |> str_subset(file.name.base) |> str_remove("\\./")
    if(length(similar.names)>1){
      latest.version <- similar.names |> str_subset("v\\d") |> str_extract("\\d$") |> as.numeric() |> max()
    } else {latest.version <- 1}
    location <- paste0(file.name.base, '---v', latest.version+1)
    folder=paste0(folder, '---v', latest.version+1)
    paste0("Simulation results for these parameters already exist. Saving results to folder ", location, ".") |> message()
  }

  dir.create(folder)
  save(output, file=paste0(folder, '/sim-data.RData'))
  env_bind(output, date=Sys.Date(),
              simulation.path=simulation.path, note=note, ec=ec)
  render('../results-template.Rmd', output_file=paste0(folder, '/report.pdf'), envir=output)
}
