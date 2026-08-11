## -----------------------------------------------------
## This script reproduces the results of
## Section B.3 of the Supplement
## -----------------------------------------------------

cat(timestamp(),'\n')

library(parallel)

## Input -----------------------------------------------------------------------

filescript <- normalizePath('Factors_V16.R', mustWork = TRUE)
source(filescript)
fileout    <- 'Section_B3.RData' # 
set.seed(212)
nwork  <- 8                # number of nodes
nrep   <- 1e3              # MC replications 
sd0    <- 1                # standard deviation of F0 (iid case)
sdE    <- 1                # standard deviation of E 
nset   <- c(20,50,100,200) # length of the series
lnset  <- length(nset)
p1set  <- c(10,20,50,100)  # rows of X 
lp1set <- length(p1set)
p2     <- 20              # columns of X
a0     <- 1               # R0,C0 ~ Uniform[-a0,a0]  + mu
a1     <- 1               # R1,C1 ~ Uniform[-a1,a1]  + mu
mu     <- 0               # mean of the above for R0,C0,R1,C1
inv.method <- 'solve'      # inversion method
phi    <- 0               # parameter for the VAR(1) Phi matrix  

## ----------------------------------------------------------------------------

hnames <-  c('hr0','hc0','hr0p','hc0p','hr1','hc1','hr1p','hc1p','hr1p2','hc1p2')
hnames <- hnames[c(1,2,5,6)]
lhnames <- length(hnames)
lhset <- 8
hset  <- matrix(c(1,1,1,1, #1
                  2,2,2,2, #2
                  1,1,2,2, #3 
                  1,1,3,3, #4
                  2,2,1,1, #5
                  3,3,1,1, #6
                  2,1,1,1, #7
                  1,1,2,1  #8 
                  ),lhset,4,byrow=TRUE)
colnames(hset) <- hnames

pnames <- c('R0','C0','R0p','C0p','R1','C1','R1p','C1p','R1p2','C1p2')
RESn <- array(NA,dim=c(nrep,lhnames,lp1set,lnset,lhset))
dimnames(RESn) <- list(1:nrep,hnames,p1set,nset,1:lhset) 

## PARAMETERS FOR PARALLEL ***************************************************
cl      <- makeCluster(nwork)
#cl      <- makeCluster(nwork,outfile='nodes_output.txt')
varlist <- c('filescript','hset','lhset','nrep','p1set','lp1set','nset','lnset','hnames', 'lhnames',
'p2','a0','a1','sd0','mu','sdE','inv.method','phi')
clusterExport(cl,varlist,envir = environment())
junk1   <- clusterEvalQ(cl, {
    source(filescript)
 })
seed    <- 123                # seed for parallel RNG
clusterSetRNGStream(cl, seed) # selects LEcuyer parallel RNG
## PARAMETERS FOR PARALLEL ***************************************************

RESn <- parSapply(cl=cl, X=1:lhset, FUN=function(ii){
    RESs <- array(NA,dim=c(nrep,lhnames,lp1set,lnset))
    dimnames(RESs) <- list(1:nrep,hnames,p1set,nset)
    hr0 <- hset[ii,1]   # rows of F0 (stationary)
    hc0 <- hset[ii,2]   # columns of F0 (stationary)
    hr1 <- hset[ii,3]   #  rows of F1 (non stationary)
    hc1 <- hset[ii,4]   # columns of F1 (non stationary)
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
          ress   <- tryCatch(FactEst(X, hr1=NULL, hc1=NULL, hr0=NULL, hc0=NULL,  hmax=hmax, inv.method = inv.method), error=function(x){list(hest=list(init=rep(NA,4),p2=rep(NA,4)))})
          RESs[h,,k,j] <- ress$hest$p2
          }
        }
    }
    return(list(ress=RESs))
}
, simplify=FALSE)



htrue <- hset 
percs <- array(NA,dim=c(lp1set,lhnames,lnset,lhset))
dimnames(percs)  <- list(p1set,hnames,nset,1:lhset)
 fnums <- mnums <- percs 
for(ii in 1:lhset){
  for(j in 1:lnset){
    for(k in 1:lp1set){
      resd <- RESn[[ii]]
      percs[k,,j,ii] <- colMeans(t(t(resd$ress[,,k,j]) == htrue[ii,]),na.rm=TRUE)*100 # percentages of correct estimation
      mnums[k,,j,ii] <- apply(resd$ress[,,k,j],MARGIN=2, FUN=median, na.rm=TRUE) # median number of factors chosen
      fnums[k,,j,ii] <- sapply(apply(resd$ress[,,k,j],MARGIN=2, FUN=tabulate,simplify=FALSE),FUN=which.max) # modal number of factors chosen
    }
  } 
}



stopCluster(cl);
#closeAllConnections();

save.image(fileout)


