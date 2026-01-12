---
title: Simulation Results
date: "2026-01-12"
output:
  pdf_document
params:
  boxplot: TRUE
  bootstrap: TRUE
  first.step: TRUE
  ci: TRUE
---








## Simulation Setup

This simulation is performed with $n=100$ and $d=40$, using the 2-d lattice as the underlying graph. $s=5$ parameters are set to be nonzero, and the beta parameter is chosen to be $\beta=0.4$. The attached results are for a 10-replication simulation. The true values of the parameter vector $\theta$ are


```
0.4472136 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4472136 0 0 0 -0.4472136 0 -0.4472136 0 -0.4472136 0 0 0 0 0 ,
```

but for brevity, our simulation only estimates the indices of $\theta$ in $\mathcal C =\{$ 1, 27, 5, 2$\}$ elements of $\theta$. Accordingly, **all statistics and visuals are indicative of performance only on the set $\mathcal C$.**







The results from our code are compared to those of [Cai, Guo, and Ma (2021)](http://www-stat.wharton.upenn.edu/~tcai/paper/Inference-GLM.pdf).





The attached results include the mean-squared error for each parameter estimate, as well as boxplots for a selection of nonzero and zero-valued parameters. In the boxplots, the green line represents the true value of the estimated parameter.

After these, I show coverage statistics for 95% symmetric confidence intervals for each of the parameters.

## Results


### Mean-squared error comparison $(\frac{1}{n.sim}\sum_{i=1}^{n.sim} \frac{1}{|\mathcal C|}\|\hat\theta_{i,\mathcal C}-\theta_{\mathcal C}\|^2)$




Table: Mean-Squared Error of Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[1]  |    5.565| 0.096|
|theta[27] |    1.494| 0.132|
|theta[5]  |    0.658| 0.051|
|theta[2]  |    7.230| 0.051|
|total     |    3.737| 0.083|



Table: Mean-Squared Error of First-Step Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[1]  |    0.280| 0.058|
|theta[27] |    0.128| 0.113|
|theta[5]  |    0.007| 0.002|
|theta[2]  |    0.000| 0.002|
|total     |    0.104| 0.044|



\newpage

### Boxplots



![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-2.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-3.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-4.pdf)<!-- --> 

\newpage

### First Step Histograms

![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-2.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-3.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n100-d40-beta0.4-s5---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-4.pdf)<!-- --> 

\newpage

### Statistics and 95% Confidence Intervals from per-Replicate Estimates


Table: Statistics for proposed Estimates

|          |    Min| Median|   Max| lower.CI.btsp| upper.CI.btsp|
|:---------|------:|------:|-----:|-------------:|-------------:|
|theta[1]  | -6.529|  0.016| 1.395|        -5.440|         1.270|
|theta[27] | -1.106|  0.843| 2.903|        -0.921|         2.747|
|theta[5]  | -0.875|  0.167| 1.575|        -0.764|         1.535|
|theta[2]  | -1.002|  0.087| 6.132|        -0.934|         6.020|


Table: Statistics for cgm Estimates

|          |    Min| Median|   Max| lower.CI.btsp| upper.CI.btsp|
|:---------|------:|------:|-----:|-------------:|-------------:|
|theta[1]  | -0.197|  0.364| 0.750|        -0.164|         0.742|
|theta[27] | -0.439|  0.311| 0.686|        -0.353|         0.683|
|theta[5]  | -0.373| -0.008| 0.403|        -0.343|         0.372|
|theta[2]  | -0.253|  0.027| 0.536|        -0.228|         0.487|

\vspace{5pt}

### Statistics for Theoretical 95% Confidence Intervals


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for proposed Estimates

|          | Estimate|    SE| lower.CI| upper.CI| cvg|
|:---------|--------:|-----:|--------:|--------:|---:|
|theta[1]  |   -0.555| 0.319|   -1.179|    0.070| 0.5|
|theta[27] |    0.868| 0.271|    0.337|    1.399| 0.4|
|theta[5]  |    0.302| 0.268|   -0.223|    0.827| 0.5|
|theta[2]  |    1.153| 0.290|    0.585|    1.721| 0.5|


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for cgm Estimates

|          | Estimate|    SE| lower.CI| upper.CI| cvg|
|:---------|--------:|-----:|--------:|--------:|---:|
|theta[1]  |    0.336| 0.170|    0.003|    0.668| 0.8|
|theta[27] |    0.265| 0.174|   -0.076|    0.606| 0.8|
|theta[5]  |   -0.013| 0.156|   -0.319|    0.294| 0.7|
|theta[2]  |    0.061| 0.163|   -0.260|    0.381| 0.9|
