# This script simulates W1 distances for one configuration from cfg.
# 
# Author: KLEY
###############################################################################

## get task id from environment variable set by SLURM:
## (if using without SLURM, set manually)
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

library(parallel)

source("../../R/AR_q.R")
source("sim_params.R")

# detach parameters from configuration
a          <- cfg[[task_id]]$a
innov      <- cfg[[task_id]]$innov
innov_name <- cfg[[task_id]]$innov_name
kappa2     <- cfg[[task_id]]$kappa2
kappa4     <- cfg[[task_id]]$kappa4
b_start    <- cfg[[task_id]]$b_start
b_end      <- cfg[[task_id]]$b_end

# Initialise L'Ecuyer's ?RngStreams?
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)

# fetch seed for this task_id's substream
s <- .Random.seed
for (i in 1:task_id) {
  s <- nextRNGStream(s)
  # send s to worker i as .Random.seed
}
.Random.seed <- s

res <- array(dim = c(R, b_end - b_start + 1, length(n_s), length(k_s) ) )


start_time <- Sys.time()
# ... then do all the computations
for (b in b_start:b_end) {
  
  perc <- 0
  for (r in 1:R) {
    
    # simulate a time series
    Y <- simAR(n = max(n_s), a = a, w = warm_up, innov = innov)
    
    for (i_n in 1:length(n_s)) {
      
      n <- n_s[i_n]

      # compute acf
      gamma <- acf(Y[1:n], lag.max = max(k_s), demean = FALSE,
                   plot = FALSE, type = "covariance")$acf[,1,1]
      
      # Code to verify what exactly acf computes:
        #> n <- length(Y)
        #> sum(Y^2)/n
        #> sum(Y[1:(n-5)]*Y[6:n])/n
      
      res[r, b - b_start + 1, i_n, ] <- gamma[which(0:max(k_s) %in% k_s)]

    }
    end_time <- Sys.time()
    if ( 100 * r / R > perc ) {
      perc <- ceiling(100 * r / R)
      cat("computed ", b, "/", B, " to ", perc, "% : ",
          difftime(end_time, start_time, units="mins"), "min\n")
    }
  }
}

# now compute the Wasserstein distance of
#  - sqrt(n) (hat gamma(k) - gamma(k)) to
#  - N(0, Sigma)

W1_gamma <- array( dim = c(b_end - b_start + 1, length(n_s), length(k_s)) )

# true gamma(k)
gamma <- acf_AR(a = a, kappa2 = kappa2)
Sigma <- Sigma_AR(a = a, kappa2 = kappa2, kappa4 = kappa4)

u <- (2*(1:R) - 1)/(2*R)
z <- qnorm( u )

start_time <- Sys.time()

for (ell in 1:(b_end - b_start + 1)) {
  for (i_n in 1:length(n_s)) {
    for (i_k in 1:length(k_s)) {
  
      X <- sqrt(n_s[i_n]) * (res[ , ell, i_n, i_k] - gamma(k_s[i_k]))
      X <- sort(X)

      s <- sqrt(Sigma(k_s[i_k], k_s[i_k]))
      W1_gamma[ell, i_n, i_k] <- mean( abs(X - s * z) )

      end_time <- Sys.time()
      cat("W1 computation: b = ", b_start + ell - 1,
          ", n = ", n_s[i_n], ", k = ", k_s[i_k], " completed (",
          difftime(end_time, start_time, units="mins"), ")\n")
    }
  }
}

save(W1_gamma,
     a, innov, innov_name,
     kappa2, kappa4, b_start, b_end,
     file = paste("W1_", task_id, ".Rdata", sep = ""))
