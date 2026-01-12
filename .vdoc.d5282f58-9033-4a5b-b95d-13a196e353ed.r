#
#
#
#
#
#
#
#
#
#
#
#
#
# render('results-template.Rmd', output_format = "pdf_document")
#
#
#
#
knitr::opts_chunk$set(echo=FALSE, comment='')

methods <- c("proposed", "cgm", "vdg")
methods.detected <- hutils::if_else(c(proposed.method, compare.to.cgm, compare.to.vdg),
                                    methods,
                                    "") |> str_subset(".")
if(!"proposed" %in% methods.detected){
  target.names <- paste0(file.param.name)
  # target.folders <- str_subset(list.files(getwd()), pattern=paste0(target.names) )
  message(list.files())
  target.folders <- grep(target.names, list.files(), value=TRUE)
  if(length(target.folders)==0) stop('No target folders found.')
  # target.folders <- list.files(pattern=fixed(target.names) )
  comparison.found <- FALSE
  for(folder in target.folders){
    candidate <- paste0(target.folders[folder], '/sim-data.RData')
    load(candidate)
    if(env_get(output, 'proposed.method')){
      results.hist <- env_get(output, 'results.hist')
      comparison.found <- TRUE
      paste0("No comparison to the proposed method was completed in the provided simulation. A comparison with the same parameters was found at ", candidate, " and will be used for comparison in this report.") |> message()
      break
    }
  }
}
#
#
#
cat(" **Note**: ", note)
#
#
#
results.env <- new_environment(list())
for(m in methods.detected){
  env_bind(results.env, !!sym(paste0("results.hist.", m)) := env_get(.GlobalEnv, paste0("results.hist.", m)) )
}
#
#
#
sqe <- matrix(-1, nrow=length(covts.record)+1, ncol=length(results.env), dimnames = list(c(paste0("theta[", covts.record, "]"), 
"total"), methods.detected))

for(m in methods.detected){
  result = env_get(results.env, paste0("results.hist.", m) )
  sqe[,m] <- result$value-matrix(theta[covts.record], byrow=TRUE, nrow=length(covts.record), ncol=n.sim)
}


for(m in methods.detected){
  for(i in 1:nrow(sqe)){
    sqe[i,as.character(m)] = (results.hist$value[i,] - theta[covts.record[i]])^2 |> mean()
    if(i==nrow(sqe)){
      sqe[i,] = mean(sqe[1:(nrow(sqe)-1),])
  }
}

}
#
#
#
#
#
#
#
cat(theta, ",")
#
#
#
#
#
#
#
#
#
cat("The results from our code are not augmented with any comparison method here.")
#
#
#
cat("The results from our code are compared to those of [Cai, Guo, and Ma (2021)](http://www-stat.wharton.upenn.edu/~tcai/paper/Inference-GLM.pdf).")
#
#
#
cat("The results from our code are compared to those of [van de Geer, Buhlmann, Ritov, and Dezeure (2014)](https://projecteuclid.org/journals/annals-of-statistics/volume-42/issue-3/On-asymptotically-optimal-confidence-regions-and-tests-for-high-dimensional/10.1214/14-AOS1221.full).")
#
#
#
cat("The results from our code are compared to those of both [Cai, Guo, and Ma (2021)](http://www-stat.wharton.upenn.edu/~tcai/paper/Inference-GLM.pdf) and [van de Geer, Buhlmann, Ritov, and Dezeure (2014)](https://projecteuclid.org/journals/annals-of-statistics/volume-42/issue-3/On-asymptotically-optimal-confidence-regions-and-tests-for-high-dimensional/10.1214/14-AOS1221.full).")
#
#
#
#
#
#
#
#
#
#
#
#
theta.mat <- matrix(0,nrow=nrow(results.hist),ncol=n.sim)
for(i in 1:nrow(results.hist)){
  theta.mat[i,] = rep(theta[results.hist$parameter.no[i]], n.sim)
}
if(params$first.step & proposed.method) tibble('MISLE (First-step) MSE' = mean((c(results.hist$first.step-theta.mat))^2), "MISLE MSE" = mean(sqe))

if(!params$first.step & proposed.method &!compare.to.cgm & !compare.to.vdg) tibble("MISLE MSE" = mean(sqe))
if(compare.to.cgm & !compare.to.vdg) tibble("MISLE MSE" = mean(sqe), "CGM MSE" = mean(sqe.cgm)) 
if(!compare.to.cgm & compare.to.vdg) tibble("MISLE MSE" = mean(sqe), "vdG MSE" = mean(sqe.vdg))
if(compare.to.cgm & compare.to.vdg) tibble("MISLE MSE" = mean(sqe), "CGM MSE" = mean(sqe.cgm), "vdG MSE" = mean(sqe.vdg))
#
#
#
#
#
cat("### Boxplots")
#
#
#
#
#
#
set.seed(1)
ind.plot <- sample(indices.on, min(s,8), replace=FALSE)
ind.plot <- c(ind.plot, sparse.ind)
for(i in ind.plot){
  data <- results.hist[results.hist$parameter.no==i,"value"]$value |> as.vector()
  # data <- results.hist[results.hist$parameter.no==i,"value"] |> as.vector()
  if(compare.to.cgm){
      data.cgm <- results.hist.cgm[results.hist.cgm$parameter.no==i,"value"]$value |> as.vector()
  # data.cgm <- results.hist.cgm[results.hist.cgm$parameter.no==i,"value"] |> as.vector()
      boxplot(list('MISLE'=data, 'CGM'=data.cgm), main=paste0("Boxplot for theta[", i, "] = ", round(theta[i], 3), "with comparison to CGM") )
  abline(h=theta[i], col='green')
  }
  if(compare.to.vdg){
      data.vdg <- results.hist.vdg[results.hist.vdg$parameter.no==i,"value"]$value |> as.vector()
  # data.vdg <- results.hist.vdg[results.hist.vdg$parameter.no==i,"value"] |> as.vector()
      boxplot(list('MISLE'=data, 'vdg'=data.vdg), main=paste0("Boxplot for theta[", i, "] = ", round(theta[i], 3), "with comparison to vdg") )
  abline(h=theta[i], col='green')
  }
}
#
#
#
#
#
#
#
cat("### First Step Histograms")
#
#
#
for(i in ind.plot){
  est <- results.hist$first.step[results.hist$parameter.no==i,]
  hist(est, freq=FALSE, main=paste0('Histogram of theta.hat[', i, ']'), xlab=NULL)
  abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
  print("Summary statistics of bootstrap replicates:")
  summary(est) |> print()
  df <- data.frame(lower=quantile(est,0.025), upper=quantile(est, 0.975))
  row.names(df) <- NULL
  print("95% CI based on bootstrap:")
  df |> print()
  if(compare.to.cgm){
    est.cgm <- results.hist.cgm$first.step[results.hist.cgm$parameter.no==i,]
    hist(est.cgm, freq=FALSE, main=paste0('Histogram of theta.hat.cgm[', i, ']'), xlab=NULL)
    abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
    print("Summary statistics of bootstrap replicates:")
    summary(est.cgm) |> print()
    df.cgm <- data.frame(lower.cgm=quantile(est.cgm,0.025), upper.cgm=quantile(est.cgm, 0.975))
    row.names(df.cgm) <- NULL
    print("95% CI based on bootstrap:")
    df.cgm |> print()
  }
}
#
#
#
#
#
#
#
for(i in ind.plot){
  est <- results.hist$value[results.hist$parameter.no==i,]
  hist(est, freq=FALSE, main=paste0('Histogram of theta.tilde[', i, ']'), xlab=NULL)
  abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
  print("Summary statistics of bootstrap replicates:")
  summary(est) |> print()
  df <- data.frame(lower=quantile(est,0.025), upper=quantile(est, 0.975))
  row.names(df) <- NULL
  print("95% CI based on bootstrap:")
  df |> print()
  if(compare.to.cgm){
    est.cgm <- results.hist.cgm$value[results.hist.cgm$parameter.no==i,]
    hist(est.cgm, freq=FALSE, main=paste0('Histogram of theta.tilde.cgm[', i, ']'), xlab=NULL)
    abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
    print("Summary statistics of bootstrap replicates:")
    summary(est.cgm) |> print()
    df.cgm <- data.frame(lower.cgm=quantile(est.cgm,0.025), upper.cgm=quantile(est.cgm, 0.975))
    row.names(df.cgm) <- NULL
    print("95% CI based on bootstrap:")
    df.cgm |> print()
  }
  if(compare.to.vdg){
    est.vdg <- results.hist.vdg$value[results.hist.vdg$parameter.no==i,]
    hist(est.vdg, freq=FALSE, main=paste0('Histogram of theta.tilde.vdg[', i, ']'), xlab=NULL)
    abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
    print("Summary statistics of bootstrap replicates:")
    summary(est.vdg) |> print()
    df.vdg <- data.frame(lower.vdg=quantile(est.vdg,0.025), upper.vdg=quantile(est.vdg, 0.975))
    row.names(df.vdg) <- NULL
    print("95% CI based on bootstrap:")
    df.vdg |> print()
  }
}
#
#
#
#
#
#
#
#
#
for(i in ind.plot){
  parameter.row <- which(results.hist$parameter.no==i)
  values <- results.hist[parameter.row, "value"]$value |> as.vector()
  ses <- results.hist[parameter.row, "se"]$se |> as.vector()
  # values <- results.hist[parameter.row, "value"] |> as.vector()
  # ses <- results.hist[parameter.row, "se"] |> as.vector()  
  intervals <- data.frame(lowers=values-1.96*ses, uppers=values+1.96*ses)
  print(paste0("Length of Confidence Intervals for theta[", i, "]"), quote=FALSE )
  print(paste0("Coverage proportion: ", mean(intervals[,1] <= theta[i] & theta[i] <= intervals[,2]) |> round(2)), quote=FALSE)
  summary(intervals$uppers-intervals$lowers) |> print(quote=FALSE)
  if(compare.to.cgm){
    values.cgm <- results.hist.cgm[parameter.row, "value"]$value |> as.vector()
    ses.cgm <- results.hist.cgm[parameter.row, "se"]$se |> as.vector()
    # values.cgm <- results.hist.cgm[parameter.row, "value"] |> as.vector()
    # ses.cgm <- results.hist.cgm[parameter.row, "se"] |> as.vector()  
    intervals.cgm <- data.frame(lowers=values.cgm-1.96*ses.cgm, uppers=values.cgm+1.96*ses.cgm)
    print(paste0("Length of Confidence Intervals for theta[", i, "]", " (CGM Method)"), quote=FALSE )  
    print(paste0("Coverage proportion: ", mean(intervals.cgm[,1] <= theta[i] & theta[i] <= intervals.cgm[,2]) |> round(2)), quote=FALSE)
    summary(intervals.cgm$uppers-intervals$lowers) |> print(quote=FALSE) 
  }
  
  if(compare.to.vdg){
    values.vdg <- results.hist.vdg[parameter.row, "value"]$value |> as.vector()
    ses.vdg <- results.hist.vdg[parameter.row, "se"]$se |> as.vector()
    # values.vdg <- results.hist.vdg[parameter.row, "value"] |> as.vector()
    # ses.vdg <- results.hist.vdg[parameter.row, "se"] |> as.vector()  
    intervals.vdg <- data.frame(lowers=values.vdg-1.96*ses.vdg, uppers=values.vdg+1.96*ses.vdg)
    print(paste0("Length of Confidence Intervals for theta[", i, "]", " (vdg Method)"), quote=FALSE )  
    print(paste0("Coverage proportion: ", mean(intervals.vdg[,1] <= theta[i] & theta[i] <= intervals.vdg[,2]) |> round(2)), quote=FALSE)
    summary(intervals.vdg$uppers-intervals$lowers) |> print(quote=FALSE) 
  }
}
#
#
#
#
