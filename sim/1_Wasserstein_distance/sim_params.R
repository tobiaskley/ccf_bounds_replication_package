# This contains the parameters with which sim/1_Wasserstein_distance/sim.R
# is run. In particular the list cfg.
# 
# Author: KLEY
###############################################################################

# how many independent replications per estimation of W1
R <- 4 * 10^6

# how many independent estimations of W1
B <- 50

warm_up <- 500

a_s <- c(0, 0.1, 0.3, 0.5, 0.7)

n_s <- c(25, 50, 75, 100, 150, 200, 250, 500, 1000, 2000)
k_s <- 0:2

innov_s <- list(function(n) {rnorm(n, 0, 1)},
                function(n) {sqrt(7/9) * rt(n, df = 9)},
                function(n) {sqrt(12/14) * rt(n, df = 14)}
                )

innov_name_s <- c("N(0,1)", "t_9", "t_14")
kappa2_s <- c(1, 1, 1)
kappa4_s <- c(0, 6/5, 6/10)

cfg <- list()

for (b in 1:B) {
    
  for (i_innov in 1:3) {
    for (a in a_s) {
        # M1
        cfg <- c(cfg, list(list(a    = a,
                                innov = innov_s[[i_innov]],
                                innov_name = innov_name_s[i_innov],
                                kappa2 = kappa2_s[i_innov],
                                kappa4 = kappa4_s[i_innov],
                                b_start = b,
                                b_end = b
                           ) ) )
  
    }
  }
}
