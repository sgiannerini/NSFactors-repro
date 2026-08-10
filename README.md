NSFactors-repro
================

This repository contains the code for replication of the results in the paper 
"Inference in High-Dimensional Matrix-Valued Time Series with Common Stochastic Trends and Multifactor Error Structure" by 
Lorenzo Trapani, Greta Goracci, Rong Chen and Simone Giannerini.


## Structure

1.  `README.md` file - This file gives a short description of the
    paper and an overview of how to carry out the analyses presented in the manuscript.
2.  `manuscript` directory - It contains 
       
       - `NSFactors_repro.Rnw`: Rnw file to reproduce all the figures and tables of both the main article and the Supplement.
       
3.  A `data` directory - This directory will generally hold the real data files 
    (or facsimile versions of them in place of confidential data) and simulated data files.
    See `data/README.md` for more details. 
4.  `code` directory - It contains

      - `Factors_V16.R`  - The library that contains the main routine `FactEst` that implements the 
         estimation methods of the paper plus auxiliary routines and `EigRatio`, that implements the
         ER criterion to estimate consistently the number of factors.
      - 
5.  An `output` directory - This directory will generally hold objects derived
    from computations, including results of simulations or real data analyses. See `output/README.md` for more details.

## Guidance on the use of reproducible environments

Submissions may include the use of reproducible environments capturing
state of a machine generating manuscript artifacts and even the
manuscript itself. Here we discuss two types of reproducible
environments and their use. Both virtual and package environments may be
put in the `code` directory.

## Language chosen

All the computations were run using `R` both under Windows and Linux. 
The typical sessionInfo is the following

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

Gentleman, Robert, and Duncan Temple Lang. “[Statistical Analyses and
Reproducible
Research](http://biostats.bepress.com/cgi/viewcontent.cgi?article=1001&context=bioconductor).”
(2004).

Gentleman, Robert. “[Reproducible research: a bioinformatics case
study](https://www.degruyter.com/document/doi/10.2202/1544-6115.1034/html).”
Statistical applications in genetics and molecular biology 4.1 (2005).

Marwick, Ben, and Bryan, Jennifer, and Attali, Dean, and Hollister,
Jeffrey W. [rrrpkg Github Page](https://github.com/ropensci/rrrpkg).
