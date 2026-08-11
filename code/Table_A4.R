## -----------------------------------------------------
## This script reproduces the results of
## Table A.4 of the Supplement
## Weak Factors case
## -----------------------------------------------------
cat(timestamp(),'\n')

library(parallel)

## Input -----------------------------------------------------------------------

filescript <- normalizePath('Factors_V16.R', mustWork = TRUE)
source(filescript)
fileout    <- 'Tabel_A4.RData' # 
set.seed(212)
nrep   <- 1e3              # MC replications 
sd0    <- 1                # standard deviation of F0 (iid case)
sdE    <- 1                # standard deviation of E 
nset   <- c(50,100,200,400,800) # length of the series
lnset  <- length(nset)
p1set  <- c(10,20,50,100)  # rows of X 
lp1set <- length(p1set)
hmax   <- 4               # max number of factors to be estimated
p2     <- 20              # columns of X
a0     <- 1               # R0,C0 ~ Uniform[-a0,a0]  + mu
a1     <- 1               # R1,C1 ~ Uniform[-a1,a1]  + mu
mu     <- 0               # mean of the above for R0,C0,R1,C1
inv.method <- 'solve'      # inversion method
phi      <- 0      # parameter for the VAR(1) Phi matrix  
alphaset <- c(0.1, 0.3, 0.5)
laset    <- length(alphaset)
nwork    <- laset                # number of nodes

## ----------------------------------------------------------------------------
make_weak_loading <- function(p, h, alpha) {
  if (h == 0L) return(matrix(0, p, 0L))
  Q <- qr.Q(qr(matrix(rnorm(p * h), p, h)))
  Q[, seq_len(h), drop = FALSE] * sqrt(p^(1 - alpha))
}
## ----------------------------------------------------------------------------

hnames <-  c('hr0','hc0','hr0p','hc0p','hr1','hc1','hr1p','hc1p','hr1p2','hc1p2')
hnames <- hnames[c(1,2,5,6)]
lhnames <- length(hnames)

hset <- c(hr0 = 1L, hc0 = 1L, hr1 = 1L, hc1 = 1L)

RESn <- array(NA, dim = c(nrep, lhnames, lp1set, lnset, laset))
dimnames(RESn) <- list(1:nrep, hnames, p1set, nset, paste0("alpha=", alphaset))

## PARAMETERS FOR PARALLEL ***************************************************
cl      <- makeCluster(nwork)
#cl      <- makeCluster(nwork,outfile='nodes_output.txt')
varlist <- c('filescript','hset','alphaset','laset','nrep','p1set','lp1set','nset','lnset',
              'hnames','lhnames','p2','a0','a1','sd0','mu','sdE','hmax',
              'inv.method','phi','make_weak_loading')
clusterExport(cl,varlist,envir = environment())
junk1   <- clusterEvalQ(cl, {
    source(filescript)
 })
seed    <- 123                # seed for parallel RNG
clusterSetRNGStream(cl, seed) # selects LEcuyer parallel RNG
## PARAMETERS FOR PARALLEL ***************************************************


 RESn <- parSapply(cl = cl, X = 1:laset, FUN = function(ia){
     alpha <- alphaset[ia]
     RESi <- array(NA, dim = c(nrep, lhnames, lp1set, lnset))
     dimnames(RESi) <- list(1:nrep, hnames, p1set, nset)
     hr0 <- unname(hset['hr0'])
     hc0 <- unname(hset['hc0'])
     hr1 <- unname(hset['hr1'])
     hc1 <- unname(hset['hc1'])
     for(j in 1:lnset){
       n <- nset[j]
#      cat('T =',n,'------------\n')
        for(k in 1:lp1set){
        p1 <- p1set[k]
        for(h in 1:nrep){
          R0 <- matrix(runif(p1*hr0,-a0,a0)+mu,p1,hr0) 
          C0 <- matrix(runif(p2*hc0,-a0,a0)+mu,p2,hc0)
          R1 <- make_weak_loading(p1, hr1, alpha)
          C1 <- make_weak_loading(p2, hc1, alpha)
          F0 <- array(NA_real_, dim = c(hr0, hc0, n))          
          F1 <- array(NA_real_, dim = c(hr1, hc1, n))
          E1 <- array(rnorm(hr1*hc1*n),dim=c(hr1,hc1,n))
          E0 <- array(rnorm(hr0*hc0*n,sd=sd0),dim=c(hr0,hc0,n)) # TODO USE rmvn
          E  <- array(rnorm(p1*p2*n,sd=sdE), dim=c(p1,p2,n))
          X  <- array(NA_real_,dim=c(p1,p2,n))
          F0[,,1] <- E0[,,1]
          F1[,,1] <- E1[,,1]
          X[,,1]  <- fit_block(R1, matrix(F1[,,1], nrow = hr1, ncol = hc1), t(C1), p1, p2) +
                     fit_block(R0, matrix(F0[,,1], nrow = hr0, ncol = hc0), t(C0), p1, p2) +
                     E[,,1]
          for(i in 2:n){
            F0[,,i] <- phi*F0[,,i-1] + E0[,,i]
            F1[,,i] <- F1[,,i-1] + E1[,,i]


            X[,,i]  <- fit_block(R1, matrix(F1[,,i], nrow = hr1, ncol = hc1), t(C1), p1, p2) +
                       fit_block(R0, matrix(F0[,,i], nrow = hr0, ncol = hc0), t(C0), p1, p2) +
                       E[,,i]
          }
          resi   <- tryCatch(FactEst(X, hr1=NULL, hc1=NULL, hr0=NULL, hc0=NULL,  hmax=hmax, inv.method = inv.method
                      ,include.zero.0=TRUE, include.zero.1=TRUE), error=function(x){list(hest=list(init=rep(NA,4),p2=rep(NA,4)))})
          RESi[h,,k,j] <- resi$hest$p2
          }
        }
    }
    return(list(resi=RESi))   
}
, simplify=FALSE)


 htrue1 <- c(hr0 = 1L, hc0 = 1L, hr1 = 1L, hc1 = 1L)
 htrue0 <- c(hr0 = 1L, hc0 = 1L, hr1 = 0L, hc1 = 0L)
 
 perci0 <- array(NA, dim = c(lp1set, lhnames, lnset, laset))
 dimnames(perci0) <- list(p1set, hnames, nset, paste0("alpha=", alphaset))
 perci1 <- perci0
 
 for(ia in 1:laset){
   for(j in 1:lnset){
     for(k in 1:lp1set){
       resd <- RESn[[ia]]
       perci0[k,,j,ia] <- colMeans(t(t(resd$resi[,,k,j]) == htrue0), na.rm = TRUE) * 100
       perci1[k,,j,ia] <- colMeans(t(t(resd$resi[,,k,j]) == htrue1), na.rm = TRUE) * 100
     }
   }
 }

stopCluster(cl);
save.image(fileout)

#closeAllConnections();
