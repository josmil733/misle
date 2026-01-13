---
title: Simulation Results
date: "2026-01-13"
output:
  pdf_document
params:
  boxplot: TRUE
  bootstrap: TRUE
  first.step: TRUE
  ci: TRUE
---








## Simulation Setup

This simulation is performed with $n=500$ and $d=10$, using the 2-d lattice as the underlying graph. $s=2$ parameters are set to be nonzero, and the beta parameter is chosen to be $\beta=0.2$. The attached results are for a 10-replication simulation. The true values of the parameter vector $\theta$ are


```
0 0.7071068 0 -0.7071068 0 0 0 0 0 0 ,
```

but for brevity, our simulation only estimates the indices of $\theta$ in $\mathcal C =\{$ 2, 4, 1, 6$\}$ elements of $\theta$. Accordingly, **all statistics and visuals are indicative of performance only on the set $\mathcal C$.**







The results from our code are compared to those of [Cai, Guo, and Ma (2021)](http://www-stat.wharton.upenn.edu/~tcai/paper/Inference-GLM.pdf).





The attached results include the mean-squared error for each parameter estimate, as well as boxplots for a selection of nonzero and zero-valued parameters. In the boxplots, the green line represents the true value of the estimated parameter.

After these, I show coverage statistics for 95% symmetric confidence intervals for each of the parameters.

## Results


### Mean-squared error comparison $(\frac{1}{n.sim}\sum_{i=1}^{n.sim} \frac{1}{|\mathcal C|}\|\hat\theta_{i,\mathcal C}-\theta_{\mathcal C}\|^2)$




Table: Mean-Squared Error of Parameter Estimates

|         | proposed|   cgm|
|:--------|--------:|-----:|
|theta[2] |    0.027| 0.618|
|theta[4] |    0.012| 0.573|
|theta[1] |    0.009| 0.015|
|theta[6] |    0.016| 0.014|
|total    |    0.016| 0.305|



Table: Mean-Squared Error of First-Step Parameter Estimates

|         | proposed|   cgm|
|:--------|--------:|-----:|
|theta[2] |    0.018| 0.244|
|theta[4] |    0.020| 0.253|
|theta[1] |    0.009| 0.000|
|theta[6] |    0.002| 0.004|
|total    |    0.012| 0.125|



\newpage

### Boxplots



![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-11-2.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-11-3.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-11-4.pdf)<!-- --> 

\newpage

### First Step Histograms

![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-13-2.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-13-3.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n500-d10-beta0.2-s2---proposed+CGM---v2\report_files/figure-latex/unnamed-chunk-13-4.pdf)<!-- --> 

\newpage

### Statistics and 95% Confidence Intervals from per-Replicate Estimates


Table: Statistics for proposed Estimates

|         |    Min| Median|    Max| lower.CI.btsp| upper.CI.btsp|
|:--------|------:|------:|------:|-------------:|-------------:|
|theta[2] |  0.351|  0.678|  0.908|         0.386|         0.885|
|theta[4] | -0.728| -0.627| -0.499|        -0.728|        -0.505|
|theta[1] | -0.134| -0.019|  0.143|        -0.130|         0.135|
|theta[6] | -0.172|  0.001|  0.286|        -0.148|         0.254|


Table: Statistics for cgm Estimates

|         |    Min| Median|   Max| lower.CI.btsp| upper.CI.btsp|
|:--------|------:|------:|-----:|-------------:|-------------:|
|theta[2] | -0.640|  0.124| 0.416|        -0.599|         0.390|
|theta[4] | -0.394| -0.082| 0.695|        -0.381|         0.628|
|theta[1] | -0.206|  0.084| 0.178|        -0.183|         0.168|
|theta[6] | -0.183| -0.001| 0.232|        -0.179|         0.196|

\vspace{5pt}

### Statistics for Theoretical 95% Confidence Intervals


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for proposed Estimates

|         | Estimate|    SE| lower.CI| upper.CI| cvg|
|:--------|--------:|-----:|--------:|--------:|---:|
|theta[2] |    0.651| 0.147|    0.363|    0.939| 0.9|
|theta[4] |   -0.635| 0.154|   -0.937|   -0.333| 1.0|
|theta[1] |   -0.001| 0.115|   -0.227|    0.225| 1.0|
|theta[6] |    0.028| 0.109|   -0.187|    0.242| 0.9|


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for cgm Estimates

|         | Estimate|    SE| lower.CI| upper.CI| cvg|
|:--------|--------:|-----:|--------:|--------:|---:|
|theta[2] |    0.002| 0.173|   -0.337|    0.341| 0.1|
|theta[4] |   -0.028| 0.170|   -0.361|    0.305| 0.0|
|theta[1] |    0.036| 0.097|   -0.153|    0.226| 0.8|
|theta[6] |   -0.008| 0.097|   -0.199|    0.182| 0.8|
