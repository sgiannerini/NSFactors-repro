## -----------------------------------------------------
## This script reproduces the results of
## Section 5 and Section B.1 of the Supplement
## -----------------------------------------------------

cat(timestamp(),'\n')

library(parallel)

## -----------------------------------------------------------------------------
## INPUT ----------------------------------------------------------------------
## -----------------------------------------------------------------------------

filescript <- normalizePath('Factors_V16.R', mustWork = TRUE)
source(filescript)
fileout <- 'Section_5_B1.RData'
set.seed(124)
nwork  <- 8                # number of nodes
nrep   <- 1e4              # MC replications 
sdE    <- 1                # standard deviation of E 
nset   <- c(20,50,100,200) # length of the series
lnset  <- length(nset)
p1set  <- c(10,20,50,100)  # rows of X
lp1set <- length(p1set)
p2     <- 20               # columns of X
mu     <- 0                # mean of the loadings R0,C0,R1,C1
inv.method <- 'solve'      # inversion method for FactEst
phi    <- 0                # parameter for the VAR(1) Phi matrix   

## -----------------------------------------------------------------------------
## ----------------------------------------------------------------------------

parnames <- c('hr0','hc0','hr1','hc1','a0','a1','sd0')
lpset    <- 8
parset <- matrix(c(1,1,1,1,1 ,1 ,1, #1
                   1,1,1,1,1 ,1 ,2, #2
                   1,1,2,2,1 ,1 ,1, #3 
                   1,1,2,2,1 ,1 ,2, #4
                   2,2,1,1,1 ,1 ,1, #5
                   2,2,1,1,1 ,1 ,2, #6                  
                   1,1,1,1,10,10,1, #7
                   1,1,1,1,10,10,2  #8
                  ),lpset,7,byrow=TRUE)
colnames(parset) <- parnames
pnames <- c('R0','C0','R0p','C0p','R1','C1','R1p','C1p','R1p2','C1p2')
lp     <- length(pnames)
RES <- array(NA,dim=c(nrep,lp,lp1set,lnset,lpset))
dimnames(RES)  <- list(1:nrep,pnames,p1set,nset,1:lpset)

## PARAMETERS FOR PARALLEL ***************************************************
cl      <- makeCluster(nwork)
#cl      <- makeCluster(nwork,outfile='nodes_output.txt')
varlist <- c('filescript','nrep','p1set','lp1set','nset','lnset','pnames',
'mu','sdE','parset','p2','lp','inv.method','phi')
clusterExport(cl,varlist,envir = environment())
junk1   <- clusterEvalQ(cl, {
    source(filescript)
 })
seed    <- 123                # seed for parallel RNG
clusterSetRNGStream(cl, seed) # selects LEcuyer parallel RNG
## PARAMETERS FOR PARALLEL ***************************************************

#for(ii in 1:lpset){
RES <- parSapply(cl=cl, X=1:lpset, FUN=function(ii){
    RESi <- array(NA,dim=c(nrep,lp,lp1set,lnset))
    dimnames(RESi) <- list(1:nrep,pnames,p1set,nset)
    hr0 <- parset[ii,1]   # rows of F0 (stationary)
    hc0 <- parset[ii,2]   # columns of F0 (stationary)
    hr1 <- parset[ii,3]   # rows of F1 (non stationary)
    hc1 <- parset[ii,4]   # columns of F1 (non stationary)
    a0  <- parset[ii,5]   # R0,C0 ~ Uniform[-a0,a0]  + mu
    a1  <- parset[ii,6]   # R1,C1 ~ Uniform[-a1,a1]  + mu
    sd0 <- parset[ii,7]   # standard deviation of F0 (iid case)
#    Phi <- diag(phi,hc0)  # VAR(1) matrix parameter
    for(j in 1:lnset){
      n <- nset[j]
#      cat('T =',n,'------------\n')
        for(k in 1:lp1set){
        p1 <- p1set[k]
        for(h in 1:nrep){
          R0 <- matrix(runif(p1*hr0,-a0,a0)+mu,p1,hr0) 
          C0 <- matrix(runif(p2*hc0,-a0,a0)+mu,p2,hc0)
          R1 <- matrix(runif(p1*hr1,-a1,a1)+mu,p1,hr1) 
          C1 <- matrix(runif(p2*hc1,-a1,a1)+mu,p2,hc1)
          F0 <- array(NA_real_, dim = c(hr0, hc0, n))          
          F1 <- array(NA_real_, dim = c(hr1, hc1, n))
          E1 <- array(rnorm(hr1*hc1*n),dim=c(hr1,hc1,n))
          E0 <- array(rnorm(hr0*hc0*n,sd=sd0),dim=c(hr0,hc0,n)) # TODO USE rmvn
          E  <- array(rnorm(p1*p2*n,sd=sdE), dim=c(p1,p2,n))
          X  <- array(NA_real_,dim=c(p1,p2,n))
          F0[,,1] <- E0[,,1]
          F1[,,1] <- E1[,,1]
          X[,,1]  <- R1 %*% F1[,,1] %*% t(C1) +  R0 %*% F0[,,1] %*% t(C0) + E[,,1]
          for(i in 2:n){
            F0[,,i] <- phi*F0[,,i-1] + E0[,,i]
            F1[,,i] <- F1[,,i-1] + E1[,,i]
            X[,,i] <- R1 %*% F1[,,i] %*% t(C1) +  R0 %*% F0[,,i] %*% t(C0) + E[,,i]
          }
          res  <- tryCatch(FactEst(X,hr1=hr1,hc1=hc1,hr0=hr0,hc0=hc0, inv.method = inv.method, fix.refined = TRUE), error=function(x){list(eigMR1=NA)})
          if(!is.na(res$eigMR1[1]))  RESi[h,,k,j] <- Dspace.Fact.n(res,R0,C0,R1,C1)
          }
        }
    }
    return(RESi)
}
, simplify='array')

mRES  <- apply(RES,MARGIN=4:2,FUN=mean) # rows = n, cols = p1
sRES  <- apply(RES,MARGIN=4:2,FUN=sd)   # rows = n, cols = p1
mcRES <- apply(RES,MARGIN=5:2,FUN=mean) # cases x  n x p1 x nparam
casenames <- paste('case',rep(1:4,each=2),rep(1:2,2),sep='.') 
dimnames(mcRES)[[1]] <- casenames

cat(timestamp(),'\n')

stopCluster(cl);
#closeAllConnections();

save(mRES,mcRES,sRES,nset,lnset,p1set,lp1set,parset,lpset,file=fileout)



