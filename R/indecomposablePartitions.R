# Returns a vector of length N+1, with the Bell numbers B_0, B_1, ..., B_N

bellNumbersUpTo <- function(N) {
  # Define B_0
  bellNumbers <- c(1)
  # Save B_n to the array at position n+1
  for (n in 0:(N-1)) {
    bellNumbers[(n+1)+1] <- sum(choose(n,0:n)*bellNumbers[1:(n+1)])
  }
  bellNumbers
}

bellNumber <- function(N) {
  bellNumbersUpTo(N)[N+1]
}

setPartitions <- function(N) {
  
  r <- rep(1,N)
  i <- N
  
  result <- c()
  
  while (prod(r != 1:N) == 0) {
    result <- c(result,r)
    if (r[i] <= max(r[1:(i-1)])) {
      r[i] <- r[i]+1
    } else {
      while (i >= 2 && (r[i] > max(r[1:(i-1)])))
      {i <- i-1}
      r[i] <- r[i]+1
      r[(i+1):N] <- rep(1,N-i)
      i <- N
    }
  }
  
  matrix(result, ncol=N, byrow=T)
}


# returns a vector (a_0, ..., a_n) with the coefficients from the binary representation
#
# sum_{i=0}^n a_i 2^i

toBinary <- function(x,n) {
  r <- x
  res <- c()
  while(x > 0) {
    res <- c(res, x %% 2)
    x <- x %/% 2
  }
  c(res,rep(0,n-length(res)))
}

verifyPartition <- function(P, J) {
  
  I <- length(J)
  N <- sum(J)
  
  indecomp <- TRUE
  
  for (n in 1:(2^I-2)) {
    M <- which(toBinary(n,I)==1)
    X <- c()
    for (m in M) {
      if (m>1) {s<-sum(J[1:m-1])+1} else {s<-1}
      e <- sum(J[1:m])
      X <- c(X,P[s:e])
    }
    Mbar <- which(toBinary(n,I)==0)
    Xbar <- c()
    for (m in Mbar) {
      if (m>1) {s<-sum(J[1:m-1])+1} else {s<-1}
      e <- sum(J[1:m])
      Xbar <- c(Xbar,P[s:e])
    }
    
    if (length(intersect(X,Xbar)) == 0) {
      indecomp <- FALSE
    }
  }
  
  indecomp
}



# this functions returns the list of all indecomposable partitions of
# the table
#
#   (1,1) ... (1,k1)
#   ...
#   (p,1) ... (p,kp)
#
# where the parameter k_s = c(k1, ..., kp) specifies the table;
# the entries k_j are positive integers.
#
# The min_size and max_size can be used to exclude partitions where any of the
# sets included fall below or above these thresholds
#
# The output is a matrix with each row representing one partition and each
# column representing one of the elements in the table, e.g. the k1+1 st column
# corresponds to element (2,1), the k1+k2+1 st column corresponds to element
# (3,1), etc.
# The entries in the table then indicate which set the element belongs to; e.g.
# if the entry in the second column equals 3, then (1,2) \in \nu_3 where
# \nu = \nu_1 \cup ... \cup \nu_R is the partition.

all_indecomposable_partitions <- function(k_s, min_size = 1,
                                          max_size = sum(k_s)) {

  n <- sum(k_s)
  B <- bellNumbersUpTo(n)[n + 1]
  P <- matrix(rep(0, (n + 1) * B), ncol = (n + 1))
  P[,1:n] <- setPartitions(n)
  
  for (i in 1:B) {
    P[i,n+1] <- verifyPartition(P[i,1:n], k_s)
    if (P[i,n+1] == 1) {
      for (v in 1:max(P[i,1:n])) {
        if (length(P[i,1:n][P[i,1:n] == v]) < min_size) {P[i,n+1] <- 0}
        if (length(P[i,1:n][P[i,1:n] == v]) > max_size) {P[i,n+1] <- 0}
      }
    }
  }
  
  Q <- matrix(P[, 1:(n + 1)][P[, n + 1] == 1], ncol = n + 1)[, 1:n]
  
  return(Q)
}

# all indecomposable partitions of the scheme
#  1 2
#  3 4
#  5 6
#  7 8
# minTwo <- all_indecomposable_partitions(c(2, 2, 2, 2), min_size = 2)

