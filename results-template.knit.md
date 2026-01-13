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

This simulation is performed with $n=200$ and $d=20$, using the 2-d lattice as the underlying graph. $s=3$ parameters are set to be nonzero, and the beta parameter is chosen to be $\beta=0.2$. The attached results are for a 10-replication simulation. The true values of the parameter vector $\theta$ are


```
0 0 0 0 0 -0.5773503 0 0 0 0 -0.5773503 0 0 -0.5773503 0 0 0 0 0 0 ,
```

but for brevity, our simulation only estimates the indices of $\theta$ in $\mathcal C =\{$ 6, 11, 4, 8$\}$ elements of $\theta$. Accordingly, **all statistics and visuals are indicative of performance only on the set $\mathcal C$.**







The results from our code are compared to those of [Cai, Guo, and Ma (2021)](http://www-stat.wharton.upenn.edu/~tcai/paper/Inference-GLM.pdf).





The attached results include the mean-squared error for each parameter estimate, as well as boxplots for a selection of nonzero and zero-valued parameters. In the boxplots, the green line represents the true value of the estimated parameter.

After these, I show coverage statistics for 95% symmetric confidence intervals for each of the parameters.

## Results


### Mean-squared error comparison $(\frac{1}{n.sim}\sum_{i=1}^{n.sim} \frac{1}{|\mathcal C|}\|\hat\theta_{i,\mathcal C}-\theta_{\mathcal C}\|^2)$




Table: Mean-Squared Error of Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[6]  |    0.108| 0.034|
|theta[11] |    0.880| 0.109|
|theta[4]  |    0.091| 0.037|
|theta[8]  |    0.132| 0.018|
|total     |    0.303| 0.049|



Table: Mean-Squared Error of First-Step Parameter Estimates

|          | proposed|   cgm|
|:---------|--------:|-----:|
|theta[6]  |    0.113| 0.050|
|theta[11] |    0.285| 0.157|
|theta[4]  |    0.022| 0.003|
|theta[8]  |    0.038| 0.006|
|total     |    0.114| 0.054|



\newpage

### Boxplots



![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-2.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-3.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-11-4.pdf)<!-- --> 

\newpage

### First Step Histograms

![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-2.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-3.pdf)<!-- --> ![](c:\Users\josmi\UFL Dropbox\Joshua Miles\Overleaf\Inference_Ising\Code\misle\n200-d20-beta0.2-s3---proposed+CGM\report_files/figure-latex/unnamed-chunk-13-4.pdf)<!-- --> 

\newpage

### Statistics and 95% Confidence Intervals from per-Replicate Estimates


Table: Statistics for proposed Estimates

|          |    Min| Median|    Max| lower.CI.btsp| upper.CI.btsp|
|:---------|------:|------:|------:|-------------:|-------------:|
|theta[6]  | -1.300| -0.782| -0.258|        -1.218|        -0.291|
|theta[11] | -0.850| -0.058|  1.053|        -0.790|         0.996|
|theta[4]  | -0.502|  0.106|  0.622|        -0.421|         0.546|
|theta[8]  | -0.727|  0.034|  0.639|        -0.644|         0.560|


Table: Statistics for cgm Estimates

|          |    Min| Median|    Max| lower.CI.btsp| upper.CI.btsp|
|:---------|------:|------:|------:|-------------:|-------------:|
|theta[6]  | -0.978| -0.586| -0.305|        -0.944|        -0.332|
|theta[11] | -0.804| -0.400|  0.108|        -0.797|         0.067|
|theta[4]  | -0.407|  0.020|  0.261|        -0.351|         0.241|
|theta[8]  | -0.093|  0.043|  0.339|        -0.083|         0.296|

\vspace{5pt}

### Statistics for Theoretical 95% Confidence Intervals


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for proposed Estimates

|          | Estimate|    SE| lower.CI| upper.CI| cvg|
|:---------|--------:|-----:|--------:|--------:|---:|
|theta[6]  |   -0.729| 0.233|   -1.186|   -0.272| 0.9|
|theta[11] |    0.112| 0.301|   -0.477|    0.701| 0.6|
|theta[4]  |    0.082| 0.231|   -0.372|    0.535| 0.9|
|theta[8]  |   -0.016| 0.218|   -0.445|    0.412| 0.8|


Table: Theoretical 95% Confidence Interval Statistics (averaged across replications) for cgm Estimates

|          | Estimate|    SE| lower.CI| upper.CI| cvg|
|:---------|--------:|-----:|--------:|--------:|---:|
|theta[6]  |   -0.599| 0.150|   -0.894|   -0.305| 0.9|
|theta[11] |   -0.388| 0.174|   -0.729|   -0.047| 0.7|
|theta[4]  |   -0.010| 0.133|   -0.271|    0.252| 0.9|
|theta[8]  |    0.066| 0.133|   -0.195|    0.327| 0.9|
