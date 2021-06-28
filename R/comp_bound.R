Q <- all_indecomposable_partitions(c(2, 2, 2, 2), min_size = 2)

## auxiliary Cpp function to speed up computing cumulants of products

Rcpp::cppFunction("
NumericMatrix comp_D(NumericVector kappa_s, IntegerMatrix Q,
                       int k, double a, int m) {

  IntegerVector dims = Q.attr(\"dim\");
  int rows_Q = dims[0]; // number of rows in Q;
  
  
  IntegerVector indices {k, 0, k, 0, 0, 0, 0, 0};
  int size_res = 2 * m + 2 * k + 1;
  NumericMatrix res = NumericMatrix(size_res, size_res);
  
  for (int u1 = 0; u1 < size_res; ++u1) {
    indices[4] = u1 - m;       // corresponds to u1 + k [-(m+k)]
    indices[5] = u1 - m - k;   // corresponds to u1     [-(m+k)]
    for (int u2 = u1; u2 < size_res; ++u2) {
      indices[6] = u2 - m;     // corresponds to u2 + k [-(m+k)]
      indices[7] = u2 - m - k; // corresponds to u2     [-(m+k)]
      double this_res = 0;
      // next four lines are to exclude cases where cumulant is 0 from
      // the costly computation
      std::vector<int> v {0, 0, 0};
      v[0] = 0; v[1] = u1; v[2] = u2;
      std::sort(v.begin(), v.end());
      if (v[1] - v[0] <= m + k or v[2] - v[1] <= m + k) {
        for (int rw = 0; rw < rows_Q; ++rw) {
          double this_fct = 1;
          IntegerVector rowQ = Q.row(rw);
          int R = max(rowQ);
          for (int r = 1; r <= R; ++r) {
            int p = 0;   // number of indices in current set
            int Min = 1000000000; // min of indices in current set
            int Max = -1000000000; // max of indices in current set
            int S = 0;   // sum of indices in current set
            
            IntegerVector::iterator iQ, iI;
            
            for (iQ = rowQ.begin(), iI = indices.begin(); iQ != rowQ.end(); ++iQ, ++iI) {
  
              if (*iQ == r) {
                p++;
                S += *iI;
                if (*iI < Min) Min = *iI;
                if (*iI > Max) Max = *iI;
              }
            }
            if (Max - Min <= m) {
              if (a == 0) {
                this_fct = 0;
              } else {
                this_fct *= kappa_s(p - 1) * pow(a, S - p * Min) *
                    (1 - pow(a, p * (m - Max + Min + 1))) / (1 - pow(a, p));
              }
            } else {
              this_fct = 0;
            }
            
          }
          this_res += this_fct;
        }
      }
      res(u1, u2) = this_res;
      res(u2, u1) = this_res;
    }
  }

  return res;
  
}")

# Function to compute the bound for given 
#   a, n, k, cumulants of epsilon_t, and
#   m = 1, ..., m_max
comp_bound_AR1 <- function(n_s, a, k, kappa_s, m_incr = 1, m_max = Inf) {

################################################################################
#  Auxiliary, private functions
################################################################################
  
  # This function computes cum(X_u1, ..., X_up) for the case when {X_t} is the
  # m-dependent approximation 
  #
  #   X_t = sum_{j=0}^m a^j eps_{t-j}
  #
  # to an AR(1) process with parameter with {eps_t} iid and cum_p(eps_t) = 1.
  #
  # Note: if cum_p(eps_t) = kappa_p, then the cumulant can be obtained by
  #       -> kappa_p * cum_AR1(u_s, a, m)
  #
  # params:
  #   u_s      vector (u_1, ..., u_p)
  #   a        autoregressive parameter in (-1,1)
  #   m        length of dependence; Inf => the AR(1) process
  #
  # returns: a vector of length n containing (X_1, ..., X_n)
  cum_AR1 <- function(u_s, a, m) {
    p <- length(u_s)
    M <- min(u_s)
    R <- max(u_s) - M
    S <- sum(u_s - M)
    if (R <= m) {
      res <- a^S * (1 - a^(p * (m - R + 1))) / (1 - a^p)
    } else {
      res <- 0
    }
    return( res )
  }

################################################################################

# For better efficiency, remove all partitions from Q that have a nu_r with
# kappa_{|nu_r|} = 0

  Q0 <- Q
  for (i in nrow(Q0):1) {
    R <- max(Q0[i,])
    rm <- FALSE
    for (r in 1:R) {
      if (kappa_s[length(which(Q0[i,] == r))] == 0) {rm <- TRUE}
    }
    if (rm) {Q0 <- Q0[-i,]}
  }

################################################################################

  # The auxiliary function here is defined for given n, a, k, kappa2, .., kappa8.
  # It computes the value of the bound for one m.
  comp_bound_aux <- function(m, n, Dm, Cm) {
    
    # Compute Sigma_k and \tilde\Sigma_k
    Sigma_k <- kappa_s[2]^2 *
      (1 + a^2 + a^(2 * abs(k)) * (1 + a^2 + 2 * k * (1 - a^2))) /
      (1 - a^2)^3 +
      kappa_s[4] * a^(2 * abs(k)) / (1 - a^2)^2
    
    # Compute \tilde\Sigma_k
    # addend <- c()
    # for (u in (-min(n - k, m)):(min(n - k, m))) {
    #   addend <- c(addend, (1 - abs(u) / (n-k)) *
    #                 (kappa_s[2]^2 * cum_AR1(c(u+k, 0), a, m) * cum_AR1(c(u-k, 0), a, m) +
    #                    kappa_s[2]^2 * cum_AR1(c(u, 0), a, m)^2 +
    #                    kappa_s[4] * cum_AR1(c(u + k, u, k, 0), a, m)) )
    # }
    # tildeSigma_k <- sum(addend) * (n - k) / n
    
    addends <- c(Cm[1])
    if (min(n - k, m) > 0) {
      addends <- c(addends,
                   2 * (1 - 1:min(n - k, m) / (n-k)) * Cm[1:min(n - k, m) + 1])
    }
    tildeSigma_k <- sum(addends) * (n - k) / n
    
    # N <- min(n - k, m)
    # tildeSigma_k <- ((n - k) / n) *
    #   (Cm[1] +
    #      2 * sum((1 - abs(1:N) / (n - k)) * Cm[2:(N + 1)]) * (N > 0) )

    # S_mn is the bound for sum tilde Q_t
    S_mn <- 0
    plus_s <- c()
    for (t in 1:(n - k)) {
      
      if (t <= 2 * (m + k) + 1 | t >= n - 3*k - 2*m) {
        
        Bt <- seq(max(1, t - 2 * (m + k)), min(t + 2 * (m + k), n - k))
        
        N <- length(Bt)
        if (N == 1) {
          M2bt <- N * Cm[1]  
        } else if (N >= 2) {
        M2bt <- N * Cm[1] +
          2 * sum( (N - (1:(N-1)) ) * Cm[2:N] )
        }
      }
      
      if (t <= (m + k) + 1 | t >= n - 2*k - m) {
        
        At <- seq(max(1, t - (m + k)), min(t + (m + k), n - k))
        
        s <- min(At) - t + (m+k) + 1
        e <- max(At) - t + (m+k) + 1
        
        M1t <- sum(Dm[s:e, s:e])

        N <- length(At)
        if (N == 1) {
          M2at <- N * Cm[1]
        } else if (N >= 2) {
          M2at <- N * Cm[1] +
                2 * sum( (N - (1:(N-1)) ) * Cm[2:N] )  * (N >= 2)
        }
        
        M3t <- 0
        for (u in At) {
            M3t <- M3t + Cm[abs(t - u) + 1]
        }
      }
      
      plus <-
        sqrt(M1t + Cm[1] * M2at + M3t^2) * sqrt(M2bt) +
        0.5 * sqrt(M1t + Cm[1] * M2at + 2 * M3t^2) * sqrt(M2at)
      plus_s <- c(plus_s, plus)
      S_mn <- S_mn + plus

    }
    
    tildeS_mn <- 2 * S_mn / (tildeSigma_k*n)^(3/2)
    
    gam_k <- kappa_s[2] * cum_AR1(c(k, 0), a, Inf)
    
    K2m <- (2 * a^(m + 1) + a^(2*m + 2)) * kappa_s[2] / (1 - a^2)
    
    return(k * abs(gam_k) / sqrt(n) +
             sqrt(2/(pi * Sigma_k)) * abs(Sigma_k - tildeSigma_k) +
             2*(n-k)*K2m/sqrt(n) + tildeS_mn)
  }

  bound_s <- list()
  min_idx <- c()
  
  # start with only m = 0
  m_s <- 0
  m <- m_s[1]
  Cm <- c()
  for (u in 0:max(4*(m+k), min(max(n_s) - k, m))) {
    Cm[u + 1] <-
      kappa_s[4]   * cum_AR1(c(k - u, -u, k, 0), a, m) +
      kappa_s[2]^2 * (cum_AR1(c(u, 0), a, m)^2 +
                        cum_AR1(c(k - u, 0), a, m) * cum_AR1(c(k + u, 0), a, m))
  }
  
  for (i_n in 1:length(n_s)) {
    min_idx[i_n] <- 1
    Dm <- comp_D(kappa_s, Q0, k, a, 0)
    bound_s[[i_n]] <- c(comp_bound_aux(m_s, n_s[i_n], Dm, Cm))
  }

  # then extend by max_incr until the best bound is not at the last m
  min_found <- rep(FALSE, length(n_s))
  while (!all(min_found) & max(m_s) <= m_max) {
    new_m_s <- (max(m_s) + 1):(max(m_s)+m_incr)
    m_s <- c(m_s, new_m_s)
    for (m in new_m_s) {
      Dm <- comp_D(kappa_s, Q0, k, a, m)
      
      Cm <- c()
      for (u in 0:max(4*(m+k), min(max(n_s) - k, m))) {
        Cm[u + 1] <-
          kappa_s[4]   * cum_AR1(c(k - u, -u, k, 0), a, m) +
          kappa_s[2]^2 * (cum_AR1(c(u, 0), a, m)^2 +
                            cum_AR1(c(k - u, 0), a, m) * cum_AR1(c(k + u, 0), a, m))
      }
      
      for (i_n in 1:length(n_s)) {
        bound_s[[i_n]] <- c(bound_s[[i_n]],
                            comp_bound_aux(m, n_s[i_n], Dm, Cm))
        min_idx[i_n] <- max(which(bound_s[[i_n]] == min(bound_s[[i_n]])))
        min_found[i_n] <- (max(m_s) > m_s[min_idx[i_n]])
      }
    }
    
    cat("computing: ", new_m_s, "\n")
  }
  
  res <- list(m_s = m_s, results = list())
  
  for (i in 1:length(min_idx)) {
    res$results[[i]] <-
      list(bound_s = bound_s[[i]],
           m_star = m_s[min_idx[i]],
           bound_star = bound_s[[i]][min_idx[i]])
  }
  return(res)
  
}
