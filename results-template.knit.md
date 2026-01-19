---
title: Simulation Results
date: "2026-01-19"
output:
  pdf_document
params:
  mse: TRUE
  mad: TRUE
  boxplot: TRUE
  bootstrap: TRUE
  first.step: TRUE
  ci: TRUE
---








## Simulation Setup

This simulation is performed with $n=200$ and $d=100$, using the 2-d lattice as the underlying graph. $s=5$ parameters are set to be nonzero, and the beta parameter is chosen to be $\beta=0.2$. The attached results are for a 10-replication simulation. The true values of the parameter vector $\theta$ are


```
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0.4472136 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0.4472136 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4472136 0 0 0 0 0 0 0 0.4472136 0 0 0 0 0 0 0 0 0.4472136 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ,
```

but for brevity, our simulation only estimates the indices of $\theta$ in $\mathcal C =\{$ 20, 51, 71, 40$\}$ elements of $\theta$. Accordingly, **all statistics and visuals are indicative of performance only on the set $\mathcal C$.**







The results from our code are compared to those of [Cai, Guo, and Ma (2021)](http://www-stat.wharton.upenn.edu/~tcai/paper/Inference-GLM.pdf).





The attached results include the mean-squared error for each parameter estimate, as well as boxplots for a selection of nonzero and zero-valued parameters. In the boxplots, the green line represents the true value of the estimated parameter.

After these, I show coverage statistics for 95% symmetric confidence intervals for each of the parameters.

## Results






```
### Mean-squared error comparison $(\frac{1}{n.sim}\sum_{i=1}^{n.sim} \frac{1}{|\mathcal C|}\|\hat\theta_{i,\mathcal C}-\theta_{\mathcal C}\|^2)$
```



Table: Mean-Squared Error of Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[20] |    0.040| 0.050|
|theta[51] |    0.071| 0.035|
|theta[71] |    0.015| 0.028|
|theta[40] |    0.020| 0.041|
|total     |    0.036| 0.038|



Table: Mean-Squared Error of First-Step Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[20] |    0.131| 0.049|
|theta[51] |    0.117| 0.035|
|theta[71] |    0.000| 0.005|
|theta[40] |    0.000| 0.003|
|total     |    0.062| 0.023|



\newpage


```
### Mean absolute deviation comparison $(\frac{1}{n.sim}\sum_{i=1}^{n.sim} \frac{1}{|\mathcal C|}\|\hat\theta_{i,\mathcal C}-\theta_{\mathcal C}\|)$
```



Table: Mean Absolute Deviation of Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[20] |    0.168| 0.180|
|theta[51] |    0.233| 0.169|
|theta[71] |    0.108| 0.126|
|theta[40] |    0.087| 0.170|
|total     |    0.149| 0.161|



Table: Mean Absolute Deviation of First-Step Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[20] |    0.352| 0.189|
|theta[51] |    0.317| 0.151|
|theta[71] |    0.006| 0.033|
|theta[40] |    0.005| 0.023|
|total     |    0.170| 0.099|



\newpage

### Boxplots



![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-12-1.pdf)<!-- --> ![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-12-2.pdf)<!-- --> ![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-12-3.pdf)<!-- --> ![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-12-4.pdf)<!-- --> 

\newpage

### First Step Histograms

![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-14-1.pdf)<!-- --> ![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-14-2.pdf)<!-- --> ![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-14-3.pdf)<!-- --> ![](C:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d100-beta0.2-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-14-4.pdf)<!-- --> 

\newpage

### Statistics and 95% Confidence Intervals from per-Replicate Estimates


Table: Statistics for proposed Estimates

|          |    Min| Median|    Max| lower.CI.btsp| upper.CI.btsp|
|:---------|------:|------:|------:|-------------:|-------------:|
|theta[20] | -0.427| -0.277| -0.084|        -0.421|        -0.097|
|theta[51] | -0.465| -0.177| -0.042|        -0.446|        -0.046|
|theta[71] | -0.248| -0.057|  0.137|        -0.221|         0.136|
|theta[40] | -0.357|  0.024|  0.053|        -0.333|         0.052|


Table: Statistics for cgm Estimates

|          |    Min| Median|    Max| lower.CI.btsp| upper.CI.btsp|
|:---------|------:|------:|------:|-------------:|-------------:|
|theta[20] | -0.989| -0.406| -0.275|        -0.919|        -0.276|
|theta[51] | -0.771| -0.429| -0.192|        -0.748|        -0.205|
|theta[71] | -0.189| -0.023|  0.417|        -0.181|         0.346|
|theta[40] | -0.332|  0.024|  0.252|        -0.332|         0.252|

\vspace{5pt}

### Statistics for Theoretical 95% Confidence Intervals


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for proposed Estimates

|          | Estimate|    SE| lower.CI| upper.CI| cvg|
|:---------|--------:|-----:|--------:|--------:|---:|
|theta[20] |   -0.279| 0.142|   -0.557|   -0.001| 0.9|
|theta[51] |   -0.218| 0.131|   -0.474|    0.038| 0.5|
|theta[71] |   -0.019| 0.131|   -0.276|    0.238| 1.0|
|theta[40] |   -0.040| 0.129|   -0.292|    0.212| 0.9|


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for cgm Estimates

|          | Estimate|    SE| lower.CI| upper.CI| cvg|
|:---------|--------:|-----:|--------:|--------:|---:|
|theta[20] |   -0.489| 0.129|   -0.743|   -0.235| 0.9|
|theta[51] |   -0.452| 0.132|   -0.711|   -0.192| 0.9|
|theta[71] |    0.013| 0.133|   -0.249|    0.274| 0.9|
|theta[40] |    0.002| 0.141|   -0.273|    0.278| 0.8|
