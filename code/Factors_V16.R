## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## Non stationary factors routines
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
##  Simone Giannerini, Greta Goracci, Lorenzo Trapani
##  May 2024
##  last modified: June 2026
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
##  LEGEND: -------------------------------------------------------------------
## ----------------------------------------------------------------------------
##  In the whole library it is assumed that
##  X is a (p1 x p2 x n) array where
##  p1: row dimension
##  p2: col dimension 
##  n : sample size  (time)
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## Convenience functions called by the main routine
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
## Helpers to allow for the zero-factor case
## ----------------------------------------------------------------------------

  zcols <- function(h) {
    if (is.null(h) || h <= 0L) integer(0) else seq_len(h)
  }
  
  safe_inv <- function(X, inv.fun, tol = NULL) {
    nr <- nrow(X); nc <- ncol(X)
    if (nr == 0L || nc == 0L) return(matrix(0, nr, nc))
    do.call(inv.fun, list(X = X, tol = tol))
  }
  
  take_eigvecs <- function(eig, h, scale = 1) {
    if (h <= 0L) return(matrix(0, nrow(eig$vectors), 0L))
    as.matrix(eig$vectors[, zcols(h), drop = FALSE] * scale)
  }
  
  orth_complement <- function(L, p, inv.fun, tol = NULL) {
    if (ncol(L) == 0L) return(diag(p))
    diag(p) - L %*% safe_inv(crossprod(L), inv.fun, tol) %*% t(L)
  }
  
  fit_block <- function(R, F, Ct, p1, p2) {
    if (ncol(R) == 0L || nrow(F) == 0L || ncol(F) == 0L || nrow(Ct) == 0L) {
      return(matrix(0, p1, p2))
    }
    R %*% F %*% Ct
  }

## ----------------------------------------------------------------------------
## Matrix computations
## ----------------------------------------------------------------------------

  accum_right_tcross <- function(X, A = NULL, transpose_x = FALSE) {
    ## Computes sum_t Z_t Z_t^T
    ## Cases:
    ##   1. Z_t = X_t
    ##   2. Z_t = X_t %*% A
    ##   3. Z_t = t(X_t)
    ##   4. Z_t = t(X_t) %*% A 
    
    dx <- dim(X)
    n  <- dx[3]
    
    out_dim <- if (transpose_x) dx[2] else dx[1]
    M <- matrix(0, nrow = out_dim, ncol = out_dim)
    
    if (!transpose_x && is.null(A)) {
      for (i in 1:n) {
        x <- X[,,i]
        M <- M + tcrossprod(x)
      }
    }
    
    if (!transpose_x && !is.null(A)) {
      for (i in 1:n) {
        x <- X[,,i]
        z <- x %*% A
        M <- M + tcrossprod(z)
      }
    }
    
    if (transpose_x && is.null(A)) {
      for (i in 1:n) {
        x <- X[,,i]
        M <- M + crossprod(x)
      }
    }
    if (transpose_x && !is.null(A)) {
      for (i in 1:n) {
        x <- X[,,i]
        z <- t(x) %*% A
        M <- M + tcrossprod(z)
      }
    }
    return(M)
  }

## -----------------------------------------------------------------------------
## Matrix inversion wrappers
## ----------------------------------------------------------------------------
  
  Chol.solve <- function(X,tol){if(is.null(tol)){tol=-1};chol2inv(chol(X,tol=tol))}                    ## Cholesky decomposition
  R.solve    <- function(X,tol){if(is.null(tol)){tol=.Machine$double.eps};solve(X,tol=tol)}            ## standard R inversion
  QR.solve   <- function(X,tol){if(is.null(tol)){tol=1e-7};qr.solve(X,tol=tol)}                        ## QR decomposition
  MP.solve   <- function(X,tol){if(is.null(tol)){tol=sqrt(.Machine$double.eps)};MASS::ginv(X,tol=tol)} ## Moore-Penrose inverse

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
## MAIN ROUTINE ---------------------------------------------------------------
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------


FactEst <- function(X, hr1=NULL, hc1=NULL, hr0=NULL, hc0=NULL, hmax=NULL, 
                 inv.method = c('solve','qr','chol','mp'), tol=NULL, forceSym=TRUE, 
                 include.zero.0 = FALSE, include.zero.1 = FALSE, fix.refined = FALSE, ...){

## Inference for factor models for matrix-valued time series, with (possibly) common 
## stochastic trend and a stationary factor structure in the error term.
## Estimates the number of factors with the Eigenvalue Ratio (ER)
## it is possible to specify an initial number of I(0) and I(1) factors, hr0, hc0, hr1, hc1. 
## If they are NULL then the ER criterion is used to estimate them. 
## This allows for iterative refinement of the estimates of the number of factors
##
## -----------------------------------------------------------------------------
## INPUT  : --------------------------------------------------------------------
## -----------------------------------------------------------------------------
## X       : matrix of time series (array with dim c(p1,p2,n))
## hr1     : initial number of rows of F1 (non stationary)    if NULL it uses the ER
## hc1     : initial number of columns of F1 (non stationary) if NULL it uses the ER
## hr0     : initial number of rows of F0 (stationary)        if NULL it uses the ER 
## hc0     : initial number of columns of F0 (stationary)     if NULL it uses the ER
## hmax    : max number of factors, if NULL defaults to min(p1,p2,n)
##           it is used for the ER   
## inv.method : matrix inversion method, can be one of the following
##              solve - standard R solve
##                 qr - uses QR decomposotion (qr.solve)
##               chol - uses Cholesky decomposotion (chol2inv) 
##                 mp - Moore-Penrose pseudo inverse (MASS::ginv)
## tol        : tolerance to be passed to inv.method, if NULL keeps the default
## forceSym   : if TRUE forces the covariance matrices to be symmetric using 
##               Matrix::forceSymmetric
## include.zero.0 :
## include.zero.1 : if TRUE includes the zero in the ER criterion
##                  for estimating the number of I(0) and I(1) factors, resp. 
## fix.refined    : if TRUE uses the user specified hr1, hc1, hr0, hc0 for hr1p, hc1p, hr0p, hc0p, hr1p2, hc1p2
##                  if FALSE (default) hr1p, hc1p, hr0p, hc0p, hr1p2, hc1p2 are always re-estimated via ER           
## ...   : further arguments passed to other routines called inside the function
##
## -----------------------------------------------------------------------------
## OUTPUT : A list with the following components: ------------------------------
## -----------------------------------------------------------------------------
##
## R1hat     : initial estimate of R1 
## C1hat     : initial estimate of C1
## eigMR1    : eigenvalues of MR1
## eigMC1    : eigenvalues of MC1
## R1phat    : projected estimate of R1
## C1phat    : projected estimate of C1
## eigMR1p   : eigenvalues of MR1p
## eigMC1p   : eigenvalues of MC1p
## R1p2hat   : iterated projected estimate of R1
## C1p2hat   : iterated projected estimate of C1
## eigMR1p2  : eigenvalues of MR1p2
## eigMC1p2  : eigenvalues of MC1p2
## R0hat     : initial estimate of R0
## C0hat     : initial estimate of C0
## eigMR0    : eigenvalues of MR0
## eigMC0    : eigenvalues of MC0
## R0phat    : final estimate of R0
## C0phat    : final estimate of C0
## eigMR0p   : eigenvalues of MR0p
## eigMC0p   : eigenvalues of MC0p
## F1p2hat   : iterated estimate of F1
## F1phat    : estimate of F1
## F0hat     : initial estimate of F0
## F0phat    : final estimate of F0
## fit       : fitted values based on F1phat and F0hat
## residuals : residuals based on F1phat and F0hat (TEMPORARILIY REMOVED)
## fitp      : fitted values based on F1phat and F0phat
## residualsp: residuals based on F1phat and F0phat (TEMPORARILIY REMOVED)
## hest      : a list with the following components
##   init    : initial number of factors  hr0  hc0  hr1  hc1 
##   p1      : refined number of factors  hr0p hc0p hr1p hc1p
##   p2      : refined number of factors  hr0p2 hc0p2 hr1p2 hc1p2
## -----------------------------------------------------------------------------
  
  inv.method   <- match.arg(inv.method)
  if(inv.method=='solve'){
    inv.fun <- 'R.solve'  
  } else
  if(inv.method=='qr'){
    inv.fun <- 'QR.solve'  
  } else
  if(inv.method=='chol'){
    inv.fun <- 'Chol.solve'  
  } else
  if(inv.method=='mp'){
    inv.fun <- 'MP.solve'  
  }

  ## Input validation ---------------------------------------------------------
 
  if (!is.array(X) || length(dim(X)) != 3L) {
    stop("'X' must be a numeric 3D array with dim = c(p1, p2, n).")
  }
  if (!is.numeric(X)) {
    stop("'X' must be numeric.")
  }
  if (any(!is.finite(X))) {
    stop("'X' contains NA, NaN or Inf values.")
  }

  dx <- dim(X)
  p1 <- dx[1]  # rows of X
  p2 <- dx[2]  # columns of X
  n  <- dx[3]  # sample size

  if (n < 2L) {
    stop("'X' must contain at least 2 time points.")
  }

check_count <- function(h, name, max_dim) {
  if (is.null(h)) return(invisible(NULL))
  if (length(h) != 1L || !is.numeric(h) || is.na(h) ||
      h != as.integer(h) || h < 0L) {
    stop(sprintf("'%s' must be a single nonnegative integer.", name))
  }
  if (h > max_dim) {
    stop(sprintf("'%s' must be <= %d.", name, max_dim))
  }
  invisible(NULL)
}

  check_count(hr0, "hr0", p1)
  check_count(hc0, "hc0", p2)
  check_count(hr1, "hr1", p1)
  check_count(hc1, "hc1", p2)

  if (is.null(hmax)) hmax <- min(p1, p2, n) # max possible number of factors
  
  if (length(hmax) != 1L || !is.numeric(hmax) || is.na(hmax) ||
      hmax != as.integer(hmax) || hmax < 2L || hmax > min(p1, p2, n)) {
    stop("'hmax' must be an integer between 2 and min(p1, p2, n).")
  }
  
#  if (xor(is.null(sc1), is.null(sr1))) {
#    stop("'sc1' and 'sr1' must be both NULL or both specified.")
#  }
  # ---------------------------------------------------------------------------

  hr0p  <- if (fix.refined) hr0 else NULL
  hc0p  <- if (fix.refined) hc0 else NULL
  hr1p  <- if (fix.refined) hr1 else NULL
  hc1p  <- if (fix.refined) hc1 else NULL
  hr1p2 <- if (fix.refined) hr1 else NULL
  hc1p2 <- if (fix.refined) hc1 else NULL
  
  # quantities for estimating the number of factors ---------------------------
  deltaR1   <- deltaC1 <- 1/n
  deltaR1p  <- 1/(sqrt(min(p1,p2))*n^(3/2)) + 1/(p2*n) + 1/n^2 + 1/(sqrt(p1*p2)*n)
  deltaC1p  <- 1/(sqrt(min(p1,p2))*n^(3/2)) + 1/(p1*n) + 1/n^2 + 1/(sqrt(p1*p2)*n)
  deltaR1p2 <- 1/(sqrt(p2)*n) + 1/n^2 +  1/(sqrt(p1)*n^(3/2))
  deltaC1p2 <- 1/(sqrt(p1)*n) + 1/n^2 +  1/(sqrt(p2)*n^(3/2))
  deltaR0   <- 1/(sqrt(p2*n)) + 1/p1 + 1/n
  deltaC0   <- 1/(sqrt(p1*n)) + 1/p2 + 1/n
  deltaR0p  <- 1/(sqrt(p2*n)) + 1/(sqrt(p1)*n) + 1/p1 + 1/n^2
  deltaC0p  <- 1/(sqrt(p1*n)) + 1/(sqrt(p2)*n) + 1/p2 + 1/n^2  
  # ------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------- --------
  ## Step 0 Preliminary estimation (flattened+projected) of the non stationary part ----
  # ------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------
      
  MR1 <- accum_right_tcross(X, A = NULL, transpose_x = FALSE) / (p1 * p2 * n^2)
  MC1 <- accum_right_tcross(X, A = NULL, transpose_x = TRUE)  / (p1 * p2 * n^2)

  if(forceSym){ 
    MR1 <- as.matrix(Matrix::forceSymmetric(MR1))
    MC1 <- as.matrix(Matrix::forceSymmetric(MC1))
  }
  eigMR1 <- eigen(MR1,symmetric=TRUE)
  eigMC1 <- eigen(MC1,symmetric=TRUE) 
  
  if(is.null(hr1)) hr1 <- EigRatio(eigMR1$values,delta=deltaR1,hmax, include.zero=include.zero.1)$index
  if(is.null(hc1)) hc1 <- EigRatio(eigMC1$values,delta=deltaC1,hmax, include.zero=include.zero.1)$index

  R1hat <- take_eigvecs(eigMR1, hr1, sqrt(p1))
  C1hat <- take_eigvecs(eigMC1, hc1, sqrt(p2))


  # ------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------
  ## Step 1 Initial Estimation of R0, C0, F0 (anti-projection based)
  # ------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------  

  C1o <- orth_complement(C1hat, p2, inv.fun, tol) # p2xp2
  R1o <- orth_complement(R1hat, p1, inv.fun, tol) # p1xp1 
  
 ## Estimation of orthogonal spaces -------------------------
 ## Stream the second moments directly instead of storing XC / XR
 
  ## NOTE: MR0 and MC0 below are called MR1o MC1o in the paper
  MR0 <- accum_right_tcross(X, C1o, transpose_x = FALSE) / (p1*p2*n)
  MC0 <- accum_right_tcross(X, R1o, transpose_x = TRUE)  / (p1*p2*n)
  
  if(forceSym){ 
    MR0 <- as.matrix(Matrix::forceSymmetric(MR0))
    MC0 <- as.matrix(Matrix::forceSymmetric(MC0))
  }
  ## Eq (9) ----------------------------
  
  eigMR0 <- eigen(MR0,symmetric=TRUE)
  eigMC0 <- eigen(MC0,symmetric=TRUE) 
  
  
  if(is.null(hr0)) hr0 <- EigRatio(eigMR0$values,delta=deltaR0, hmax, include.zero=include.zero.0)$index
  if(is.null(hc0)) hc0 <- EigRatio(eigMC0$values,delta=deltaC0, hmax, include.zero=include.zero.0)$index
  

  R0hat <- take_eigvecs(eigMR0, hr0, sqrt(p1))
  C0hat <- take_eigvecs(eigMC0, hc0, sqrt(p2))
  
  tR0hat <- t(R0hat)
  tC0hat <- t(C0hat)

  # L0 is part of the vectorised estimator of F0 ------------------------------

  A1 <- tC0hat%*%C1o%*%C0hat
  B1 <- tR0hat%*%R1o%*%R0hat
  
  A1i <- safe_inv(A1, inv.fun, tol)
  B1i <- safe_inv(B1, inv.fun, tol)
  
  K1 <- kronecker(A1i, B1i, make.dimnames = FALSE)  
  K2 <- kronecker((tC0hat%*%C1o), (tR0hat%*%R1o), make.dimnames = FALSE)
  L0 <- K1 %*% K2  ## first two factors of Eq.(10) 

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ## Step 2 Projection-based Estimation of R1, C1, F1 
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  F0hat <- array(0, dim = c(hr0, hc0, n))
  MR1p  <- matrix(0, nrow = p1, ncol = p1)
  MC1p  <- matrix(0, nrow = p2, ncol = p2)

  ## pass 1: estimate F0hat and accumulate MR1p / MC1p without storing Xp
 for(i in 1:n){
  x <- X[,,i]
  if (hr0 == 0L || hc0 == 0L) {
      F0i <- matrix(0, hr0, hc0)
    } else {
      F0i <- matrix(L0 %*% matrix(x), hr0, hc0)
    }
    F0hat[,,i] <- F0i

    xf <- x - fit_block(R0hat, F0i, tC0hat, p1, p2)

    xfc <- xf %*% C1hat
    xfr <- t(xf) %*% R1hat
    MR1p <- MR1p + tcrossprod(xfc)
    MC1p <- MC1p + tcrossprod(xfr)
  }
 
  MR1p <- MR1p / (p1  * p2^2* n^2)
  MC1p <- MC1p / (p1^2* p2  * n^2)
  
  if(forceSym){ 
    MR1p <- as.matrix(Matrix::forceSymmetric(MR1p))
    MC1p <- as.matrix(Matrix::forceSymmetric(MC1p))
  }
  
  ## Eq.(12,13) ----------------------------
  eigMR1p <- eigen(MR1p,symmetric=TRUE)
  eigMC1p <- eigen(MC1p,symmetric=TRUE) 
  
  if(is.null(hr1p)) hr1p <- EigRatio(eigMR1p$values,delta=deltaR1p, hmax, include.zero=include.zero.1)$index
  if(is.null(hc1p)) hc1p <- EigRatio(eigMC1p$values,delta=deltaC1p, hmax, include.zero=include.zero.1)$index
  

  R1phat <- take_eigvecs(eigMR1p, hr1p, sqrt(p1)) # refined estimation of R1
  C1phat <- take_eigvecs(eigMC1p, hc1p, sqrt(p2)) # refined estimation of C1
  
  tR1phat <- t(R1phat)
  tC1phat <- t(C1phat)
  
  ## pass 2: estimate F1phat by recomputing the filtered slice on the fly
  if (hr1p == 0L || hc1p == 0L) {
    F1phat <- array(0, dim = c(hr1p, hc1p, n))
  } else {
    F1phat <- array(0, dim = c(hr1p, hc1p, n))
    for (i in 1:n) {
      x <- X[,,i] 
  
      if (hr0 == 0L || hc0 == 0L) {
        xf <- x
      } else {
        xf <- x - R0hat %*% matrix(F0hat[,,i], nrow = hr0, ncol = hc0) %*% tC0hat
      }
      F1phat[,,i] <- tR1phat %*% xf %*% C1phat
    }
    F1phat <- F1phat / (p1 * p2)
  }
 
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ## Step 3 iterated anti-projection/projection-based Estimation of R0, C0, F0 
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
 ## FINAL Estimation of R0, C0, F0 -------------------------------
  
 
  C1po <- orth_complement(C1phat, p2, inv.fun, tol) # p2xp2
  R1po <- orth_complement(R1phat, p1, inv.fun, tol) # p1xp1 
   
  ## FINAL estimation of the I(0) loadings R0 C0 ----
  ## Stream the second moments directly instead of storing XCp / XRp
  
  ## NOTE: MR0p and MC0p below are called MR1op MC1op in the paper
  MR0p <- accum_right_tcross(X, C1po%*%C0hat, transpose_x = FALSE) / (p1  *p2^2*n)
  MC0p <- accum_right_tcross(X, R1po%*%R0hat, transpose_x = TRUE)  / (p1^2*p2  *n)
  
 if(forceSym){ 
    MR0p <- as.matrix(Matrix::forceSymmetric(MR0p))
    MC0p <- as.matrix(Matrix::forceSymmetric(MC0p))
  }
  # Add Eq. to the paper
  eigMR0p <- eigen(MR0p,symmetric=TRUE)
  eigMC0p <- eigen(MC0p,symmetric=TRUE) 
  
  ### PUNTO DA CHIEDERE A LORENZO hr0p va stimato o è =hr0?
  if(is.null(hr0p)) hr0p <- EigRatio(eigMR0p$values,delta=deltaR0p, hmax, include.zero=include.zero.0)$index
  if(is.null(hc0p)) hc0p <- EigRatio(eigMC0p$values,delta=deltaC0p, hmax, include.zero=include.zero.0)$index
  

  R0phat <- take_eigvecs(eigMR0p, hr0p, sqrt(p1))
  C0phat <- take_eigvecs(eigMC0p, hc0p, sqrt(p2))
  
  tR0phat <- t(R0phat)
  tC0phat <- t(C0phat)
  
  ## Iterated Estimation of F0 -----------------------
  # L0p is the time-independent part of the vectorised estimator of F0p
#  if((hr0==hr0p)&(hc0==hc0p)){
    A1p <- tC0phat%*%C1po%*%C0phat
    B1p <- tR0phat%*%R1po%*%R0phat
    K2p <- kronecker((tC0phat%*%C1po), (tR0phat%*%R1po), make.dimnames = FALSE)
#  }else{
#    A1p <- tC0phat%*%C1pos%*%C0phat
#    B1p <- tR0phat%*%R1pos%*%R0phat
#    K2p <- kronecker((tC0phat%*%C1pos), (tR0phat%*%R1pos), make.dimnames = FALSE)
#  }
  
  A1pi <- safe_inv(A1p, inv.fun, tol)
  B1pi <- safe_inv(B1p, inv.fun, tol)
  
  K1p <- kronecker(A1pi, B1pi, make.dimnames = FALSE)  
 
  L0p <- K1p %*% K2p 
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ## Step 4 iterated projection-based estimation of R1, C1, F1 
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  MR1p2 <- matrix(0, nrow = p1, ncol = p1)
  MC1p2 <- matrix(0, nrow = p2, ncol = p2)

## pass 1: estimate F0phat and accumulate MR1p2 / MC1p2
  F0phat <- array(0, dim = c(hr0p, hc0p, n))

  for (i in 1:n) {
    x <- X[,,i]
    if (hr0p == 0L || hc0p == 0L) {
      xf <- x
    } else {
      F0i <- matrix(L0p %*% matrix(x), hr0p, hc0p) ##  Eq.(15)
      F0phat[,,i] <- F0i
      xf <- x - R0phat %*% F0i %*% tC0phat  # filtered Xt
    }
    xfc <- xf %*% C1phat
    xfr <- t(xf) %*% R1phat
    MR1p2 <- MR1p2 + tcrossprod(xfc)
    MC1p2 <- MC1p2 + tcrossprod(xfr)
  }  

  MR1p2 <- MR1p2 / (p1  * p2^2* n^2)
  MC1p2 <- MC1p2 / (p1^2* p2  * n^2)
  
  if(forceSym){ 
    MR1p2 <- as.matrix(Matrix::forceSymmetric(MR1p2))
    MC1p2 <- as.matrix(Matrix::forceSymmetric(MC1p2))
  }
  
  ## Eq.(?) ----------------------------
  eigMR1p2 <- eigen(MR1p2,symmetric=TRUE)
  eigMC1p2 <- eigen(MC1p2,symmetric=TRUE) 

#  TO BE DISCUSSED  
  if(is.null(hr1p2)) hr1p2 <- EigRatio(eigMR1p2$values,delta=deltaR1p2, hmax, include.zero=include.zero.1)$index
  if(is.null(hc1p2)) hc1p2 <- EigRatio(eigMC1p2$values,delta=deltaC1p2, hmax, include.zero=include.zero.1)$index
  
 
  R1p2hat <- take_eigvecs(eigMR1p2, hr1p2, sqrt(p1)) # final estimation of R1 
  C1p2hat <- take_eigvecs(eigMC1p2, hc1p2, sqrt(p2)) # final estimation of C1
 
  tR1p2hat <- t(R1p2hat)
  tC1p2hat <- t(C1p2hat)
  
  ## pass 2: estimate F1p2hat by recomputing the filtered slice on the fly

  if (hr1p2 == 0L || hc1p2 == 0L) {
    F1p2hat <- array(0, dim = c(hr1p2, hc1p2, n))
  } else {
    F1p2hat <- array(0, dim = c(hr1p2, hc1p2, n))
  for (i in 1:n) {
      x <- X[,,i]
      if (hr0p == 0L || hc0p == 0L) {
        xf <- x
      } else {
xf <- x - R0phat %*% matrix(F0phat[,,i], nrow = hr0p, ncol = hc0p) %*% tC0phat
      }
      F1p2hat[,,i] <- tR1p2hat %*% xf %*% C1p2hat
    }
    F1p2hat <- F1p2hat / (p1 * p2)   ## Eq(14)
  }  
 
 
  ## Predicted values -----------------------
  Xhat <- Xphat <- Xp2hat  <- array(NA,dim=c(p1,p2,n))
  for(i in 1:n){
    Xhat[,,i] <- fit_block(R1phat,  matrix(F1phat[,,i],  nrow = hr1p,  ncol = hc1p),  tC1phat,  p1, p2) +
             fit_block(R0hat,   matrix(F0hat[,,i],   nrow = hr0,   ncol = hc0),   tC0hat,   p1, p2)

    Xphat[,,i] <- fit_block(R1phat,  matrix(F1phat[,,i],  nrow = hr1p,  ncol = hc1p),  tC1phat,  p1, p2) +
              fit_block(R0phat,  matrix(F0phat[,,i],  nrow = hr0p,  ncol = hc0p),  tC0phat,  p1, p2)

    Xp2hat[,,i] <- fit_block(R1p2hat, matrix(F1p2hat[,,i], nrow = hr1p2, ncol = hc1p2), tC1p2hat, p1, p2) +
               fit_block(R0phat,  matrix(F0phat[,,i],  nrow = hr0p,  ncol = hc0p),  tC0phat,  p1, p2)
  }
   
  res <- list(
  R1hat  =R1hat,  C1hat  =C1hat,  eigMR1  =eigMR1$values,  eigMC1  =eigMC1$values, 
  R1phat =R1phat, C1phat =C1phat, eigMR1p =eigMR1p$values, eigMC1p =eigMC1p$values,
  R1p2hat=R1p2hat,C1p2hat=C1p2hat,eigMR1p2=eigMR1p2$values,eigMC1p2=eigMC1p2$values,
  R0hat  =R0hat , C0hat  =C0hat , eigMR0  =eigMR0$values , eigMC0  =eigMC0$values,
  R0phat =R0phat, C0phat =C0phat, eigMR0p =eigMR0p$values, eigMC0p =eigMC0p$values,
  F1phat =F1phat, F1p2hat =F1p2hat, F0hat=F0hat, F0phat=F0phat, 
  fit = Xhat, fitp = Xphat, fitp2 = Xp2hat, 
  #residuals = (X - Xhat), residualsp = X - Xphat,
  hest = list( init=c(hr0=unname(hr0), hc0=unname(hc0), hr1=unname(hr1), hc1=unname(hc1)),
               p1=c(hr0p=unname(hr0p), hc0p=unname(hc0p),hr1p=unname(hr1p), hc1p=unname(hc1p)),  
               p2=c(hr0p=unname(hr0p), hc0p=unname(hc0p),hr1p2=unname(hr1p2), hc1p2=unname(hc1p2))))  
  return(res)
} 

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------



## ----------------------------------------------------------------------------
## Other functions
## ----------------------------------------------------------------------------

DSpace <- function(A, B) {
  qa <- NCOL(A)
  qb <- NCOL(B)

  if (qa == 0L && qb == 0L) return(0)
  if (qa == 0L || qb == 0L) return(1)

  m  <- max(qa, qb)
  An <- svd(A)$u
  Bn <- svd(B)$u
  sqrt(1 - sum(diag(An %*% t(An) %*% Bn %*% t(Bn))) / m)
}
## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
  
Dspace.Fact <- function(fit,R0,C0,R1,C1){
  # Given a fit and the true parameters, it computes the distances 
  #  between true and estimated spaces
  return(c(
    R0=DSpace(R0,fit$R0hat),
    C0=DSpace(C0,fit$C0hat),
    R1=DSpace(R1,fit$R1hat),
    C1=DSpace(C1,fit$C1hat),
    R1p=DSpace(R1,fit$R1phat),
    C1p=DSpace(C1,fit$C1phat)))
}

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------

Dspace.Fact.n <- function(fit,R0,C0,R1,C1){
  # Given a fit and the true parameters, it computes the distances 
  #  between true and estimated spaces
  return(c(
    R0=DSpace(R0,fit$R0hat),
    C0=DSpace(C0,fit$C0hat),
    R0p=DSpace(R0,fit$R0phat),
    C0p=DSpace(C0,fit$C0phat),
    R1=DSpace(R1,fit$R1hat),
    C1=DSpace(C1,fit$C1hat),
    R1p=DSpace(R1,fit$R1phat),
    C1p=DSpace(C1,fit$C1phat),
    R1p2=DSpace(R1,fit$R1p2hat),
    C1p2=DSpace(C1,fit$C1p2hat)))
}

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------

EigRatio <- function(eig, delta, hmax, include.zero=TRUE){
# estimation of the number of factors using the 
# Eigenvalue Ratio
  k    <- sum(eig)
  eig  <- eig[1:hmax]
  eig0 <- eig[1:(hmax-1)]
  if(include.zero){
    omega <- sqrt(delta) # mock lambda_0
    eig0  <- c(omega,eig0)
    eig1  <- eig + k*delta
    res   <- which.max(eig0/eig1)-1
  }else{
    eig1 <- eig[2:hmax] + k*delta
    res  <- which.max(eig0/eig1)
  }
  r.max <- max(eig0/eig1)
    return(list(index=res,ratio=r.max))
}

## ----------------------------------------------------------------------------
## ----------------------------------------------------------------------------
