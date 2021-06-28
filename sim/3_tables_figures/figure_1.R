# This script returns the plots shown in Figure 1.
# 
# Author: KLEY
###############################################################################


source("../../R/indecomposablePartitions.R")
source("../../R/comp_bound.R")

## Case 1: eps_t ~ N(0,1)
kappa_s_1 <- c(0, 1, 0, 0, 0, 0, 0, 0)


## Case 2: eps_t ~ nu^{-1/2} (nu - 2)^{1/2} t_nu
nu <- 14
kappa_s_2 <- c(0, 1, 0, 6 / (nu - 4), 0, 240 / ((nu - 4) * (nu - 6)), 0,
               5040 * (5*nu - 22) / ((nu - 4)^2 * (nu - 6) * (nu - 8)))
nu <- 9
kappa_s_3 <- c(0, 1, 0, 6 / (nu - 4), 0, 240 / ((nu - 4) * (nu - 6)), 0,
               5040 * (5*nu - 22) / ((nu - 4)^2 * (nu - 6) * (nu - 8)))


# Fig 1

k <- 0
a <- 0.7

n_s <- c(25, 500, 2000)
b <- comp_bound_AR1(n_s, a, k, kappa_s_1, m_incr = 30, m_max = 30)

for (i in 1:3) {
  pdf(paste("bound_fct_N_a07_n", n_s[i], ".pdf", sep = ""), width = 6, height = 4)
    plot(b$m_s, b$results[[i]]$bound_s, type = "l",
         ylim = c(0, max(b$results[[i]]$bound_s)),
         xlab = expression(m), ylab = "value of the bound",
         main = paste("n =", n_s[i], "B_n(m*) =", round(b$results[[i]]$bound_star, 3)))
    abline(v = b$results[[i]]$m_star, lty = 2)
    abline(h = b$results[[i]]$bound_star, col = "gray")
    points(b$m_s, b$results[[i]]$bound_s)
  dev.off()
}


# Fig 2

k <- 0
a <- 0.7

n_s <- c(25, 500, 2000)
b <- comp_bound_AR1(n_s, a, k, kappa_s_3, m_incr = 30, m_max = 30)

for (i in 1:3) {  
  pdf(paste("bound_fct_t9_a07_n", n_s[i], ".pdf", sep = ""), width = 6, height = 4)
  plot(b$m_s, b$results[[i]]$bound_s, type = "l",
       ylim = c(0, max(b$results[[i]]$bound_s)),
       xlab = expression(m), ylab = "value of the bound",
       main = paste("n =", n_s[i], "B_n(m*) =", round(b$results[[i]]$bound_star, 3)))
  abline(v = b$results[[i]]$m_star, lty = 2)
  abline(h = b$results[[i]]$bound_star, col = "gray")
  points(b$m_s, b$results[[i]]$bound_s)
  dev.off()
}



# Fig 3

k <- 0
a <- 0.7

n_s <- c(25, 500, 2000)
b <- comp_bound_AR1(n_s, a, k, kappa_s_2, m_incr = 30, m_max = 30)
  
for (i in 1:3) {
  pdf(paste("bound_fct_t14_a07_n", n_s[i], ".pdf", sep = ""), width = 6, height = 4)
    plot(b$m_s, b$results[[i]]$bound_s, type = "l",
         ylim = c(0, max(b$results[[i]]$bound_s)),
         xlab = expression(m), ylab = "value of the bound",
         main = paste("n =", n_s[i], "B_n(m*) =", round(b$results[[i]]$bound_star, 3)))
    abline(v = b$results[[i]]$m_star, lty = 2)
    abline(h = b$results[[i]]$bound_star, col = "gray")
    points(b$m_s, b$results[[i]]$bound_s)
  dev.off()
}
