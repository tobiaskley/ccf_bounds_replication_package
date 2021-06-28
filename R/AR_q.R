# This script provides functions to compute the constants in the bounds for
# AR(q) processes with iid innovations
# 
# Author: KLEY
###############################################################################

# This function simulates a stretch of length n from an AR(q) process
#
# Precisely the following autoregressive scheme is used:
#
# X[t] = a[1] X[t-1] + ... + a[q]*X[t-q] + eps[t], t = -w+1, ..., n
# X[t] = eps_init[-t-w+1], t = -w-q+1, ..., -w
#
# params:
#   n        length of the stretch to be returned
#   a        vector (a_1, ..., a_q) with the coefficients of
#            the AR(q) process
#   w        length of the warm-up/burn-in period; default 100
#   innov    a function(n) that returns n innovations
#
# returns: a vector of length n containing (X_1, ..., X_n)

simAR <- function (n, a, w = 100,
                   innov = rnorm) {

  #   eps      vector with n+w innovations used for the autoregression
  #   eps_init vector with q innovations used to initialise
  eps      <- innov(w+n)
  eps_init <- innov(length(a))
  
  Y <- filter(eps, a, method = "recursive", init = eps_init)
  
  return(Y[(w+1):(w+n)])
  
}

# This function returns the transfer function A(lambda) for an AR(q) process
#
# params:
#   a      vector (a_1, ..., a_q) with the coefficients of the AR(q) process
#
# returns: a function with one argument (i.e., lambda)
trFct_AR <- function(a) {
  q <- length(a)
  A <- function(lambda) {
    1 - sum( a * exp(-1i * lambda * 1:q))
  }
  return(Vectorize(A))
}

# This function returns the power spectrum f2(lambda) for an AR(q) process
#
# params:
#   a      vector (a_1, ..., a_q) with the coefficients of the AR(q) process
#   kappa2 variance of eps_t
#
# returns: a function with one argument (i.e., lambda)
pwrSpec_AR <- function(a, kappa2) {
  A <- trFct_AR(a)
  f2 <- function(lambda) {
    abs(1/A(lambda))^2 * kappa2 / (2*pi)
  }
  return(Vectorize(f2))
}

# This function returns the tri spectrum f4(lambda1,lambda2, lambda3) for an AR(q) process
#
# params:
#   a      vector (a_1, ..., a_q) with the coefficients of the AR(q) process
#   kappa4 4th cumulant of eps_t
#
# returns: a function with three arguments (i.e., lambda1, lambda2, lambda3)
triSpec_AR <- function(a, kappa4 = 0) {
  A <- trFct_AR(a)
  f4 <- function(lambda1, lambda2, lambda3) {
    1/(A(lambda1) * A(lambda2) * A(lambda3) * 
      A(-lambda1-lambda2-lambda3)) * kappa4 / (2*pi)^3
  }
  return(Vectorize(f4))
}


# This function returns the acf for an AR(q) process
#
# params:
#   a      vector (a_1, ..., a_q) with the coefficients of the AR(q) process
#   kappa2 variance of eps_t
#
# returns: a list with the four constants as named elements
acf_AR <- function(a, kappa2 = 1) {
  
  f2 <- pwrSpec_AR(a, kappa2)

  gamma <- function (k) {
    aux <- function(alpha) {
      cos( -k * alpha) * f2(alpha)
    }
    return(2 * integrate( aux, 0, pi )$value)
  }
  
  return(Vectorize(gamma))
  
}

# This function returns the 4th order cumulant function for an AR(q) process
#
# params:
#   a      vector (a_1, ..., a_q) with the coefficients of the AR(q) process
#   kappa4 4th cumulant of eps_t
#
# returns: a list with the four constants as named elements
c4_AR <- function(a, kappa4) {
  
  f4 <- triSpec_AR(a, kappa4)
  
  c4 <- function (k1, k2, k3) {
    # Re part of Third thing
    aux_1_Re <- function(lambda1) {
      aux_2_Re <- function(lambda2) {
        aux_3_Re <- function(lambda3) {
          val <- Re( exp(1i * (k1 * lambda1 + k2 * lambda2 + k3 * lambda3) ) * 
                       f4(lambda1, lambda2, lambda3) )
          return(val)
        }
        return(integrate( Vectorize(aux_3_Re), 0, 2*pi )$value)
      }
      return(integrate( Vectorize(aux_2_Re), 0, 2*pi )$value)
    }
    return(integrate( Vectorize(aux_1_Re), 0, 2*pi )$value)
  }

  return(Vectorize(c4))
  
}


# This function returns the asymptotic covariance Sigma(k1,k2) of the acf
# when it is computed from an AR(q) process
#
# params:
#   a      vector (a_1, ..., a_q) with the coefficients of the AR(q) process
#   kappa2 variance of eps_t
#   kappa4 4th cumulant of eps_t
#
# returns: the asymptotic covariance Sigma(k1,k2)
Sigma_AR <- function(a, kappa2, kappa4) {
  
  f2 <- pwrSpec_AR(a, kappa2)
  f4 <- triSpec_AR(a, kappa4)
  
  Sigma <- function (k1, k2) {
    
    # Re part of first thing
    aux_Re <- function(alpha) {
      Re( exp(-1i * (k1 + k2) * alpha)) * f2(alpha) * f2(-alpha)
    }
    S1 <- integrate( aux_Re, 0, 2*pi )$value
    
    # Re part of second thing
    aux_Re <- function(alpha) {
      Re( exp(1i * (k1 - k2) * alpha)) * f2(alpha) * f2(-alpha)
    }
    S2 <- integrate( aux_Re, 0, 2*pi )$value
    
    # Re part of Third thing
    aux_outer_Re <- function(beta) {
      aux_inner_Re <- function(alpha) {
        Re( exp(1i * (k1 * alpha - k2 * beta) ) * f4(alpha, -alpha, -beta) )
      }
      integrate( aux_inner_Re, 0, 2*pi )$value
    }
    S3 <- integrate( Vectorize(aux_outer_Re), 0, 2*pi )$value
    
    return(2 * pi * (S1 + S2 + S3))
  }
  return(Vectorize(Sigma))
}
