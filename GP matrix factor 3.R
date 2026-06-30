###############################################################
## Gaussian Process Matrix Factor Model
###############################################################
set.seed(1)
rm(list=ls())
###############################################################################
# Read Excel workbook
# Each sheet = One observation matrix
###############################################################################

library(readxl)
library(doParallel)
library(foreach)

file <- "BrainRegion_CommonRegions.xlsx"

# Sheet names
sheet_names <- excel_sheets(file)

cat("Number of observations =", length(sheet_names), "\n\n")


###############################################################################
# Read all sheets
###############################################################################

Xlist <- vector("list", length(sheet_names))

for(i in seq_along(sheet_names)){
  
  cat("Reading", sheet_names[i], "... ")
  
  df <- read_excel(file,
                   sheet = sheet_names[i])
  
  df <- as.data.frame(df)
  
  # First column contains row names (Words)
  rownames(df) <- df[[1]]
  
  # Remove first column
  df <- df[, -1, drop = FALSE]
  
  # Convert to numeric matrix
  X <- as.matrix(df)
  
  storage.mode(X) <- "double"
  
  Xlist[[i]] <- X
  
  cat(dim(X)[1], "x", dim(X)[2], "\n")
  
}

###############################################################################
# Check dimensions
###############################################################################

dims <- t(sapply(Xlist, dim))

colnames(dims) <- c("Rows", "Columns")

print(dims)

###############################################################################
# Verify all matrices have same dimension
###############################################################################

if(length(unique(dims[,1])) != 1)
  stop("Row dimensions differ.")

if(length(unique(dims[,2])) != 1)
  stop("Column dimensions differ.")

cat("\nAll matrices have identical dimensions.\n")

###############################################################################
# Names
###############################################################################

row_names <- rownames(Xlist[[1]])
col_names <- colnames(Xlist[[1]])

###############################################################################
# Number of observations
###############################################################################

N <- length(Xlist)
IN <- diag(N)
###############################################################################
# Matrix dimensions
###############################################################################

p <- nrow(Xlist[[1]])
q <- ncol(Xlist[[1]])
###############################################################
## Dimensions
###############################################################
p
q
r <- 2
s <- 2
N

###############################################################
## Precompute observed vectors
###############################################################

Xarray <- array(0,c(p,q,N))

for(k in 1:N)
  Xarray[,,k] <- Xlist[[k]]

###############################################################################
## Common quantities for likelihood and gradients
###############################################################################

gp_common_quantities <- function(
    R,
    C,
    Farray,
    phi,
    tau,
    sigma2,
    Xarray
){
  
  ###############################################################
  ## Dimensions
  ###############################################################
  
  p <- dim(Xarray)[1]
  q <- dim(Xarray)[2]
  N <- dim(Xarray)[3]
  
  IN <- diag(N)
  
  ###############################################################
  ## Tarray = R F_k C'
  ###############################################################
  
  Tarray <- array(0, c(p, q, N))
  
  for(k in 1:N)
    Tarray[,,k] <- R %*% Farray[,,k] %*% t(C)
  
  ###############################################################
  ## Storage
  ###############################################################
  
  DeltaArray <- array(0, c(N, N, p, q))
  Uarray     <- array(0, c(N, N, p, q))
  SigmaArray <- array(0, c(N, N, p, q))
  KinvArray  <- array(0, c(N, N, p, q))
  AlphaArray <- array(0, c(N, p, q))
  M0Array    <- array(0, c(N, N, p, q))
  
  loglik <- 0
  
  ###############################################################
  ## Loop over (i,j)
  ###############################################################
  
  for(i in 1:p){
    
    for(j in 1:q){
      
      #######################################################
      ## Mean trajectory
      #######################################################
      
      mu <- Tarray[i,j,]
      
      #######################################################
      ## Pairwise differences
      #######################################################
      
      Delta <- outer(mu, mu, "-")
      
      U <- Delta^2
      
      #######################################################
      ## Covariance
      #######################################################
      
      Sigma <- phi * exp(-0.5 * tau * U)
      
      #######################################################
      ## Covariance matrix
      #######################################################
      
      K <- Sigma + sigma2 * IN
      
      #######################################################
      ## Cholesky
      #######################################################
      
      L <- chol(K)
      
      #######################################################
      ## K^{-1}
      #######################################################
      
      Kinv <- chol2inv(L)
      
      #######################################################
      ## alpha = K^{-1}x
      #######################################################
      
      x <- Xarray[i,j,]
      
      alpha <- backsolve(
        L,
        forwardsolve(t(L), x)
      )
      
      #######################################################
      ## M0
      #######################################################
      
      M0 <- tcrossprod(alpha) - Kinv
      
      #######################################################
      ## Log-likelihood
      #######################################################
      
      logdet <- 2 * sum(log(diag(L)))
      
      loglik <- loglik -
        0.5 * (
          N * log(2*pi) +
            logdet +
            crossprod(x, alpha)
        )
      
      #######################################################
      ## Save
      #######################################################
      
      DeltaArray[,,i,j] <- Delta
      Uarray[,,i,j]     <- U
      SigmaArray[,,i,j] <- Sigma
      KinvArray[,,i,j]  <- Kinv
      AlphaArray[,i,j]  <- alpha
      M0Array[,,i,j]    <- M0
      
    }
  }
  
  ###############################################################
  ## Return
  ###############################################################
  
  list(
    loglik     = as.numeric(loglik),
    Tarray     = Tarray,
    DeltaArray = DeltaArray,
    Uarray     = Uarray,
    SigmaArray = SigmaArray,
    KinvArray  = KinvArray,
    AlphaArray = AlphaArray,
    M0Array    = M0Array,
    IN         = IN
  )
}

###############################################################
## Initial Parameters
###############################################################

R <- qr.Q(qr(matrix(rnorm(p*r), p, r)))
C <- qr.Q(qr(matrix(rnorm(q*s), q, s)))

# R <- diag(p)[, 1:r, drop = FALSE]
# C <- diag(q)[, 1:s, drop = FALSE]

## Initial latent factors
Farray <- array(0, dim = c(r, s, N))

for(k in 1:N){
  Farray[,,k] <- t(R) %*% Xarray[,,k] %*% C
}

phi <- 1
tau <- 1
sigma2 <- 0.1

###############################################################################
## Gradient with respect to (phi, tau, sigma2)
###############################################################################

gradient_theta <- function(
    R,
    C,
    Farray,
    phi,
    tau,
    sigma2,
    Xarray
){
  
  ###############################################################
  ## Common quantities
  ###############################################################
  
  cache <- gp_common_quantities(
    R       = R,
    C       = C,
    Farray  = Farray,
    phi     = phi,
    tau     = tau,
    sigma2  = sigma2,
    Xarray  = Xarray
  )
  
  ###############################################################
  ## Gradient accumulators
  ###############################################################
  
  del_phi    <- 0
  del_tau    <- 0
  del_sigma2 <- 0
  
  ###############################################################
  ## Identity matrix
  ###############################################################
  
  IN <- cache$IN
  
  ###############################################################
  ## Loop over (i,j)
  ###############################################################
  
  for(i in 1:p){
    
    for(j in 1:q){
      
      #######################################################
      ## Cached quantities
      #######################################################
      
      Sigma <- cache$SigmaArray[,,i,j]
      U     <- cache$Uarray[,,i,j]
      M0    <- cache$M0Array[,,i,j]
      
      #######################################################
      ## Derivatives of K
      #######################################################
      
      dK_phi <- Sigma / phi
      
      dK_tau <- -0.5 * Sigma * U
      
      dK_sigma2 <- IN
      
      #######################################################
      ## Accumulate gradients
      #######################################################
      
      del_phi <-
        del_phi +
        0.5 * sum(M0 * t(dK_phi))
      
      del_tau <-
        del_tau +
        0.5 * sum(M0 * t(dK_tau))
      
      del_sigma2 <-
        del_sigma2 +
        0.5 * sum(M0 * t(dK_sigma2))
      
    }
    
  }
  
  ###############################################################
  ## Return gradient
  ###############################################################
  
  c(
    del_phi,
    del_tau,
    del_sigma2
  )
  
}

###############################################################################
## Gradient with respect to R
###############################################################################

gradient_R <- function(
    R,
    C,
    Farray,
    phi,
    tau,
    sigma2,
    Xarray
){
  
  ###########################################################################
  ## Common quantities
  ###########################################################################
  
  cache <- gp_common_quantities(
    R       = R,
    C       = C,
    Farray  = Farray,
    phi     = phi,
    tau     = tau,
    sigma2  = sigma2,
    Xarray  = Xarray
  )
  
  DeltaArray <- cache$DeltaArray
  SigmaArray <- cache$SigmaArray
  M0Array    <- cache$M0Array
  
  ###########################################################################
  ## Gradient
  ###########################################################################
  
  grad <- matrix(0, p, r)
  
  ###########################################################################
  ## Precompute F^(k) C[j,]
  ##
  ## FCarray[k,,j] = F^(k) %*% C[j,]
  ###########################################################################
  
  FCarray <- array(0, c(N, r, q))
  
  for(j in 1:q){
    
    cj <- C[j,]
    
    for(k in 1:N){
      
      FCarray[k,,j] <- Farray[,,k] %*% cj
      
    }
    
  }
  
  ###########################################################################
  ## Precompute pairwise differences
  ##
  ## DeltaFCarray[,,n,j]
  ###########################################################################
  
  DeltaFCarray <- array(0, c(N, N, r, q))
  
  for(j in 1:q){
    
    FC <- FCarray[,,j]
    
    for(n in 1:r){
      
      DeltaFCarray[,,n,j] <-
        outer(
          FC[,n],
          FC[,n],
          "-"
        )
      
    }
    
  }
  
  ###########################################################################
  ## Compute gradient
  ###########################################################################
  
  for(m in 1:p){
    
    ## only i = m contributes
    
    for(j in 1:q){
      
      ###################################################################
      ## Common part (independent of n)
      ###################################################################
      
      Base <-
        M0Array[,,m,j] *
        t(DeltaArray[,,m,j]) *
        t(SigmaArray[,,m,j])
      
      ###################################################################
      ## Loop over latent dimension
      ###################################################################
      
      for(n in 1:r){
        
        grad[m,n] <-
          grad[m,n] -
          0.5 * tau *
          sum(
            Base *
              t(DeltaFCarray[,,n,j])
          )
        
      }
      
    }
    
  }
  
  grad
  
}
###############################################################################
## Gradient with respect to C
###############################################################################

gradient_C <- function(
    C,
    R,
    Farray,
    phi,
    tau,
    sigma2,
    Xarray
){
  
  ###########################################################################
  ## Common quantities
  ###########################################################################
  
  cache <- gp_common_quantities(
    R       = R,
    C       = C,
    Farray  = Farray,
    phi     = phi,
    tau     = tau,
    sigma2  = sigma2,
    Xarray  = Xarray
  )
  
  DeltaArray <- cache$DeltaArray
  SigmaArray <- cache$SigmaArray
  M0Array    <- cache$M0Array

  ###########################################################################
  ## Gradient
  ###########################################################################
  
  grad <- matrix(0, q, s)
  
  ###########################################################################
  ## Precompute R[i,] %*% F^(k)
  ##
  ## RFarray[k,,i] = R[i,] %*% F^(k)
  ###########################################################################
  
  RFarray <- array(0, c(N, s, p))
  
  for(i in 1:p){
    
    ri <- R[i,]
    
    for(k in 1:N){
      
      RFarray[k,,i] <- ri %*% Farray[,,k]
      
    }
    
  }
  
  ###########################################################################
  ## Pairwise differences
  ##
  ## DeltaRFarray[,,n,i]
  ###########################################################################
  
  DeltaRFarray <- array(0, c(N, N, s, p))
  
  for(i in 1:p){
    
    RF <- RFarray[,,i]
    
    for(n in 1:s){
      
      DeltaRFarray[,,n,i] <-
        outer(
          RF[,n],
          RF[,n],
          "-"
        )
      
    }
    
  }
  
  ###########################################################################
  ## Compute gradient
  ###########################################################################
  
  for(m in 1:q){
    
    ## only j = m contributes
    
    for(i in 1:p){
      
      #####################################################################
      ## Common part
      #####################################################################
      
      Base <-
        M0Array[,,i,m] *
        t(DeltaArray[,,i,m]) *
        t(SigmaArray[,,i,m])
      
      #####################################################################
      ## Latent dimension
      #####################################################################
      
      for(n in 1:s){
        
        grad[m,n] <-
          grad[m,n] -
          0.5 * tau *
          sum(
            Base *
              t(DeltaRFarray[,,n,i])
          )
        
      }
      
    }
    
  }
  
  grad
  
}

###############################################################################
## Gradient with respect to F
###############################################################################

gradient_F <- function(
    Farray,
    R,
    C,
    phi,
    tau,
    sigma2,
    Xarray
){
  
  ###########################################################################
  ## Dimensions
  ###########################################################################
  
  p <- nrow(R)
  q <- nrow(C)
  
  r <- ncol(R)
  s <- ncol(C)
  
  N <- dim(Xarray)[3]
  
  ###########################################################################
  ## Common quantities
  ###########################################################################
  
  cache <- gp_common_quantities(
    R       = R,
    C       = C,
    Farray  = Farray,
    phi     = phi,
    tau     = tau,
    sigma2  = sigma2,
    Xarray  = Xarray
  )
  
  DeltaArray <- cache$DeltaArray
  SigmaArray <- cache$SigmaArray
  M0Array    <- cache$M0Array
  
  ###########################################################################
  ## Gradient
  ###########################################################################
  
  grad <- array(0, c(r, s, N))
  
  ###########################################################################
  ## Loop over observed matrix entries
  ###########################################################################
  
  for(i in 1:p){
    
    Ri <- R[i,]
    
    for(j in 1:q){
      
      Cj <- C[j,]
      
      #######################################################################
      ## Outer product:
      ## RC[m,n] = R[i,m] * C[j,n]
      #######################################################################
      
      RC <- tcrossprod(Ri, Cj)
      
      #######################################################################
      ## Common N x N matrix
      #######################################################################
      
      Base <- M0Array[,,i,j] * SigmaArray[,,i,j]
      
      #######################################################################
      ## Loop over latent observations
      #######################################################################
      
      for(k in 1:N){
        
        for(l in 1:N){
          
          if(l == k)
            next
          
          coeff <-
            -tau *
            Base[k,l] *
            DeltaArray[k,l,i,j]
          
          grad[,,k] <-
            grad[,,k] + coeff * RC
          
        }
        
      }
      
    }
    
  }
  
  grad
  
}
###############################################################################
## Block Coordinate L-BFGS
###############################################################################

max_iter <- 100
tol <- 0.1

loglik_old <- -Inf
###############################################################
## History table
###############################################################

history <- data.frame(
  Iteration = integer(),
  LogLik    = numeric(),
  phi        = numeric(),
  tau        = numeric(),
  sigma2     = numeric(),
  Rmax       = numeric(),
  Cmax       = numeric(),
  Fmax       = numeric(),
  Fnorm      = numeric()
)
for(iter in 1:max_iter){
  
  cat("\n========================================\n")
  cat("Outer Iteration :", iter, "\n")
  cat("========================================\n")
  
  ###########################################################################
  ## STEP 1 : Update (phi, tau, sigma2)
  ###########################################################################
  
  theta0 <- c(phi, tau, sigma2)
  cat("Starting theta optimization...\n")
  flush.console()
  fit_theta <- optim(
    
    par = theta0,
    
    fn = function(theta){
      
      -gp_common_quantities(
        R       = R,
        C       = C,
        Farray  = Farray,
        phi     = theta[1],
        tau     = theta[2],
        sigma2  = theta[3],
        Xarray  = Xarray
      )$loglik
      
    },
    
    gr = function(theta){
      
      -gradient_theta(
        R       = R,
        C       = C,
        Farray  = Farray,
        phi     = theta[1],
        tau     = theta[2],
        sigma2  = theta[3],
        Xarray  = Xarray
      )
      
    },
    
    method = "L-BFGS-B",   
    control = list(
      trace  = 1,
      REPORT = 1,
      maxit  = 10
    ),
    lower = c(1e-6,1e-6,1e-6)
    
  )
  
  phi    <- fit_theta$par[1]
  tau    <- fit_theta$par[2]
  sigma2 <- fit_theta$par[3]
  cat("Finished theta optimization\n")
  flush.console()
  ###########################################################################
  ## STEP 2 : Update R
  ###########################################################################
  cat("Starting R optimization...\n")
  flush.console()
  fit_R <- optim(
    
    par = as.vector(R),
    
    fn = function(rvec){
      
      Rnew <- matrix(rvec,p,r)
      
      -gp_common_quantities(
        R       = Rnew,
        C       = C,
        Farray  = Farray,
        phi     = phi,
        tau     = tau,
        sigma2  = sigma2,
        Xarray  = Xarray
      )$loglik
    },
    
    gr = function(rvec){
      
      Rnew <- matrix(rvec,p,r)
      
      as.vector(
        
        -gradient_R(
          R       = Rnew,
          C       = C,
          Farray  = Farray,
          phi     = phi,
          tau     = tau,
          sigma2  = sigma2,
          Xarray  = Xarray
        )
        
      )
      
    },
    
    method = "L-BFGS-B",  
    control = list(
      trace  = 1,
      REPORT = 1,
      maxit  = 10
    )
    
  )
  
  ###############################################################
  ## Orthogonalize R and absorb transformation into F
  ###############################################################
  
  Rtemp <- matrix(fit_R$par, p, r)
  
  QR <- qr(Rtemp)
  
  Q <- qr.Q(QR)
  
  A <- qr.R(QR)
  
  R <- Q
  
  for(k in 1:N)
    Farray[,,k] <- A %*% Farray[,,k]
  ###############################################################
  ## Sign convention for R
  ###############################################################
  
  for(j in 1:r){
    
    ind <- which.max(abs(R[,j]))
    
    if(R[ind,j] < 0){
      
      R[,j] <- -R[,j]
      
      for(k in 1:N)
        Farray[j,,k] <- -Farray[j,,k]
    }
  }
  cat("Finished R optimization\n")
  flush.console()
  ###########################################################################
  ## STEP 3 : Update C
  ###########################################################################
  cat("Starting C optimization...\n")
  flush.console()
  fit_C <- optim(
    
    par = as.vector(C),
    
    fn = function(cvec){
      
      Cnew <- matrix(cvec,q,s)
      
      -gp_common_quantities(
        R       = R,
        C       = Cnew,
        Farray  = Farray,
        phi     = phi,
        tau     = tau,
        sigma2  = sigma2,
        Xarray  = Xarray
      )$loglik
      
    },
    
    gr = function(cvec){
      
      Cnew <- matrix(cvec,q,s)
      
      as.vector(
        
        -gradient_C(
          C       = Cnew,
          R       = R,
          Farray  = Farray,
          phi     = phi,
          tau     = tau,
          sigma2  = sigma2,
          Xarray  = Xarray
        )
        
      )
      
    },
    
    method = "L-BFGS-B", 
    control = list(
      trace  = 1,
      REPORT = 1,
      maxit  = 10
    )
  )
  
  ###############################################################
  ## Orthogonalize C and absorb transformation into F
  ###############################################################
  
  Ctemp <- matrix(fit_C$par, q, s)
  
  QR <- qr(Ctemp)
  
  Q <- qr.Q(QR)
  
  B <- qr.R(QR)
  
  C <- Q
  
  for(k in 1:N)
    Farray[,,k] <- Farray[,,k] %*% t(B)
  ###############################################################
  ## Sign convention for C
  ###############################################################
  
  for(j in 1:s){
    
    ind <- which.max(abs(C[,j]))
    
    if(C[ind,j] < 0){
      
      C[,j] <- -C[,j]
      
      for(k in 1:N)
        Farray[,j,k] <- -Farray[,j,k]
    }
  }
  cat("Finished C optimization\n")
  flush.console()
  ###########################################################################
  ## STEP 4 : Update Farray
  ###########################################################################
  cat("Starting F optimization...\n")
  flush.console()
  fit_F <- optim(
    
    par = as.vector(Farray),
    
    fn = function(fvec){
      
      Farr <- array(fvec,c(r,s,N))
      
      -gp_common_quantities(
        R       = R,
        C       = C,
        Farray  = Farr,
        phi     = phi,
        tau     = tau,
        sigma2  = sigma2,
        Xarray  = Xarray
      )$loglik
    },
    
    gr = function(fvec){
      
      Farr <- array(fvec,c(r,s,N))
      
      as.vector(
        
        -gradient_F(
          Farray  = Farr,
          R       = R,
          C       = C,
          phi     = phi,
          tau     = tau,
          sigma2  = sigma2,
          Xarray  = Xarray
        )
        
      )
    },
    
    method = "L-BFGS-B", 
    control = list(
      trace  = 1,
      REPORT = 1,
      maxit  = 10
    )
  )
  
  Farray <- array(fit_F$par,c(r,s,N))
  cat("Finished F optimization\n")
  flush.console()
  ###########################################################################
  ## Current log-likelihood
  ###########################################################################
  
  loglik_new <- gp_common_quantities(
    R       = R,
    C       = C,
    Farray  = Farray,
    phi     = phi,
    tau     = tau,
    sigma2  = sigma2,
    Xarray  = Xarray
  )$loglik
  
  cat("Log-likelihood :", loglik_new, "\n")
  ###############################################################
  ## Save history
  ###############################################################
  
  history <- rbind(
    history,
    data.frame(
      Iteration = iter,
      LogLik    = loglik_new,
      phi        = phi,
      tau        = tau,
      sigma2     = sigma2,
      Rmax       = max(abs(R)),
      Cmax       = max(abs(C)),
      Fmax       = max(abs(Farray),
      Fnorm      = sqrt(sum(Farray^2)))
    )
  )
  print(history)
  ###########################################################################
  ## Convergence
  ###########################################################################
  
  if(abs(loglik_new - loglik_old) < tol){
    
    cat("Converged.\n")
    
    break
    
  }
  
  loglik_old <- loglik_new
  
}
cat("\n=====================================\n")
cat("Final History\n")
cat("=====================================\n")
print(history)