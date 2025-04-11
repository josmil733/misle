#' @export report

report <- function(simulation.path=NULL,lattice=TRUE, n=400, d=40, p.edge=NULL, s=5, beta=0.2, seed=0, n.burnin=30000, n.sim=100, keep.every=5, verbose=TRUE,
            n.lambda=20, eps = .00001, tau=0.8, sample.split=TRUE,
            compare.to.cgm=TRUE,
            note=NULL, ec=FALSE
            ){

# report(lattice=TRUE, n=20, d=4, p.edge=0.5, s=2, beta=0.05, seed=0, n.burnin=2000, n.sim=2, keep.every=5, verbose=TRUE, n.lambda=20, eps=1e-5, tau=0.8, sample.split=TRUE, compare.to.cgm = TRUE, note="this simulation is made with beta=0 to ensure that the results match those of CGM.")

# simulate(lattice=TRUE, n=20, d=4, p.edge=0.5, s=2, beta=0.05, seed=0, n.burnin=2000, n.sim=2, keep.every=5, verbose=TRUE, n.lambda=20, eps=1e-5, tau=0.8, sample.split=TRUE, compare.to.cgm = TRUE)

 # perform the simulations with varying n (=400, ... 1000), covariate dimensions (take d= n/10), sparsity level (s=2,5,10). Change beta around some as well. Repeat experiments B=50 times and report what proportion of confidence intervals cover the true parameters. Write a description along with the plots of all the results you received

require(rlang)
require(stringr)
source("C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/misle-estimation.R")

if(!is.null(simulation.path)){
  load(simulation.path)
} else {
output <- simulate(lattice, n, d, p.edge, s, beta, seed, n.burnin, n.sim, keep.every, verbose,
           n.lambda, eps, tau, sample.split, compare.to.cgm)
}

save(output, file=paste0("C:/Users/josmi/UFL Dropbox/Joshua Miles/Overleaf/Inference_Ising/Code/hold---", "beta", env_get(output, 'beta'), "-s", env_get(output,'s'), ".RData") ) #temporary measure

# sim.code <- runif(1, 999, 9999) |> ceiling() #to identify data and report
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
  env_bind(output,
           params=list(date=Sys.Date(),
              note=note,
              compare.to.cgm=compare.to.cgm)
  )
  # render('results-template.Rmd', output_file=paste0(folder, '/report.pdf')
  render('../results-template.Rmd', output_file=paste0(folder, '/report.pdf'), envir=output)
}
