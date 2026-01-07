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
library(knitr)
opts_chunk$set(echo=FALSE, comment='')

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
results.env <- new_environment(list())
for(m in methods.detected){
  env_bind(results.env, !!sym(paste0("results.hist.", m)) := env_get(.GlobalEnv, paste0("results.hist.", m)) )
}
#
#
#
sqe <- matrix(NA, nrow=length(covts.record)+1, ncol=length(results.env), dimnames = list(c(paste0("theta[", covts.record, "]"), 
"total"), methods.detected))

sqe.first.step = sqe

for(m in methods.detected){
  result = env_get(results.env, paste0("results.hist.", m) )
  sqe[-nrow(sqe),m] <- (result$value-matrix(theta[covts.record], nrow=length(covts.record), ncol=n.sim))^2 |> rowMeans()
  sqe[nrow(sqe),m] <- mean(sqe[-nrow(sqe),m])
  if(m != 'vdg'){
    sqe.first.step[-nrow(sqe),m] <- (result$first.step-matrix(theta[covts.record], nrow=length(covts.record), ncol=n.sim))^2 |> rowMeans()
    sqe.first.step[nrow(sqe),m] <- mean(sqe.first.step[-nrow(sqe),m])
  }
}

kable(sqe, digits=3, caption="Mean-Squared Error of Parameter Estimates")
kable(sqe.first.step, digits=3, caption="Mean-Squared Error of First-Step Parameter Estimates")
#
#
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
par(mfrow=c(1,length(methods.detected)))
for(i in covts.record){
  for(m in methods.detected){
    result = env_get(results.env, paste0("results.hist.", m) )
    hist(result$value, freq = FALSE, main=paste0("Histogram of ", m, " estimates for theta[", i, "]=", theta[i]), cex.main=0.6, xlab=NULL)
    abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
  }
}
#
#
#
#
#
cat("### First Step Histograms")
#
#
#
par(mfrow=c(1,length(methods.detected)))
for(i in covts.record){
  for(m in methods.detected){
    result = env_get(results.env, paste0("results.hist.", m) )
    hist(result$first.step, freq = FALSE, main=paste0("Histogram of ", m, " first-step estimates for theta[", i, "]=", theta[i]), cex.main=0.5, xlab=NULL)
    abline(v=theta[i], col='green', lty='dashed', lwd=2.5)
  }
}
#
#
#
#
#
#
#
out = list()
for(m in methods.detected){
  result = env_get(results.env, paste0('results.hist.', m))
  env_bind(.GlobalEnv, !!sym(paste0("df.out.", m)) := data.frame())
}

for(m in methods.detected){
    result = env_get(results.env, paste0('results.hist.', m))
    result.df = env_get(.GlobalEnv, paste0('df.out.', m))
for(i in seq_along(covts.record)){
    result.augment = data.frame(Min = min(result$value[i,]), Median = median(result$value[i,]), Max = max(result$value[i,]), lower.CI.btsp=quantile(result$value[i,],0.025), upper.CI.btsp=quantile(result$value[i,],0.975))
    if(i==1) result.df = result.augment else result.df = rbind(result.df, result.augment)
  }
  rownames(result.df) = paste0("theta[", covts.record, "]")
  kable(result.df, digits=3, caption=paste0("Bootstrap Statistics for ", m, " Estimates") )
  # env_bind(.GlobalEnv, !!sym(paste0("df.out.", m)) := result.df)
  # kable(env_get(.GlobalEnv, paste0('df.out.', m)), digits=3, caption=paste0("Bootstrap Statistics for ", m, " Estimates") )
}
#
#
#
#
#
#
#
for(m in methods.detected){
  result = env_get(results.env, paste0('results.hist.', m))
  env_bind(.GlobalEnv, !!sym(paste0("df.out.", m)) := data.frame())
}

for(m in methods.detected){
    result = env_get(results.env, paste0('results.hist.', m))
    result.df = env_get(.GlobalEnv, paste0('df.out.', m))
for(i in seq_along(covts.record)){
    result.augment = data.frame(Estimate = mean(result$value[i,]), SE = mean(result$se[i,]), lower.CI=mean(result$ci[i,,'lower']), upper.CI=mean(result$ci[i,,'upper']))
    if(i==1) result.df = result.augment else result.df = rbind(result.df, result.augment)
  }
  rownames(result.df) = paste0("theta[", covts.record, "]")
  env_bind(.GlobalEnv, !!sym(paste0("df.out.", m)) := result.df)
}
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
