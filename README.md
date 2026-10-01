NSFactors-repro
================

This repository contains the code for replication of the results in the paper 
"Inference in High-Dimensional Matrix-Valued Time Series with Common Stochastic Trends and Multifactor Error Structure" by 
Lorenzo Trapani, Greta Goracci, Rong Chen and Simone Giannerini. All the figures and tables of both the main article and the Supplement
can be reproduced by knitting the file `NSFactors_repro.Rnw`. See below for details.

We develop a factor model for high-dimensional matrix-valued time series, where both common stochastic 
trends and common stationary factors are present. We study the estimation of loading spaces, of the common stochastic trends and factors. 
Standard projection-based techniques fail to improve the rates of convergence compared to a flattened estimation technique 
which does not take into account the matrix nature of the data. Hence, we develop an iterative procedure based on anti-projections. 
We show that this estimator succeeds in refining the rates of convergence of the initial, flattened estimator. 
Monte Carlo evidence shows that our estimators offer a marked improvement compared to the flattened techniques. 
We validate our theory through an application to a set of multinational macroeconomic indices.

## Structure

1.  `README.md` file - This file gives a short description of the
    paper and an overview of how to carry out the analyses presented in the manuscript.
2.  `manuscript` directory - It contains 
       
       - `NSFactors_repro.Rnw`: Rnw file to reproduce all the figures and tables of both the main article and the Supplement.
          It uses the `RData` contained in the `output` directory, which can be produced by sourcing
          the R files contained in the `code` directory. Please take a look at the `setup` chunk to learn about the 
          configuration used to produce the figures. In particular, notice the use of `pdfcrop` and `dev='tikz'`. If you want to speed up compilation
          you can use `dev='pdf'`, instead.
       
3.  `data` directory - This directory contains the real data files
      
      - `OECDdata.csv`: 7 macroeconomic indicators for 10 OECD countries, Australia, Austria, Canada, Finland, France, Germany, New Zealand, Norway, Sweden and UK.
         Data are recorded on a quarterly basis between 2002:Q1 until 2019:Q4.
  
4.  `code` directory - It contains 

      - `Factors_V16.R`  - The workhorse library that includes 
        * `FactEst` - the main routine  that implements the estimation methods of the paper 
        * `EigRatio` - Implements the ER criterion to estimate consistently the number of factors.
        * Auxiliary routines
    
    Also, it contains the following scripts that can be sourced directly to produce the `RData` files
    included in the `output` directory:
    
      - `Section_4_B1.R`
      - `Section_B2.R`
      - `Section_B3.R`
      - `Table_A1.R`
      - `Table_A2.R`
      - `Table_A3.R`
      - `Table_A4.R`
      
5.  `output` directory - It contains the outputs of the corresponding `R` scripts contained in the `code` directory.
    The outputs are used inside `NSFactors_repro.Rnw` to produce figures and tables.

      - `Section_4_B1.RData`
      - `Section_B2.RData`
      - `Section_B3.RData`
      - `Table_A1.RData`
      - `Table_A2.RData`
      - `Table_A3.RData`
      - `Table_A4.RData`


## Language details and session info

All the computations were run using `R 4.6.0` linked against `OpenBLAS version 0.3.33` in a multicore environment, 
both under Windows and Linux. The typical sessionInfo under Windows is the following

```
R version 4.6.0 (2026-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.0
  OpenBLAS version 0.3.33 
locale:
[1] LC_COLLATE=Italian_Italy.utf8  LC_CTYPE=Italian_Italy.utf8   
[3] LC_MONETARY=Italian_Italy.utf8 LC_NUMERIC=C                  
[5] LC_TIME=Italian_Italy.utf8    

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] kableExtra_1.4.0 knitr_1.51      

loaded via a namespace (and not attached):
 [1] svglite_2.2.2      cli_3.6.6          rlang_1.2.0        xfun_0.57         
 [5] stringi_1.8.7      otel_0.2.0         textshaping_1.0.5  glue_1.8.1        
 [9] htmltools_0.5.9    scales_1.4.0       rmarkdown_2.31     evaluate_1.0.5    
[13] fastmap_1.2.0      lifecycle_1.0.5    stringr_1.6.0      compiler_4.6.0    
[17] RColorBrewer_1.1-3 rstudioapi_0.18.0  systemfonts_1.3.2  farver_2.1.2      
[21] digest_0.6.39      viridisLite_0.4.3  R6_2.6.1           magrittr_2.0.5    
[25] tools_4.6.0        xml2_1.5.2        

```


## References

Chen, R., Giannerini, S., Goracci, G., & Trapani, L. (2025). 
Inference in matrix-valued time series with common stochastic trends and multifactor error structure. 
arXiv preprint [arXiv:2501.01925](https://arxiv.org/abs/2501.01925).
