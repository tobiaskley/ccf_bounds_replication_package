# This contains the parameters with which sim/2_bound/sim.R is run.
# In particular the list cfg.
# 
# Author: KLEY
###############################################################################

## Distributions for eps_t are

## Case 1: eps_t ~ N(0,1)
kappa_s_1 <- c(0, 1, 0, 0, 0, 0, 0, 0)


## Case 2: eps_t ~ nu^{-1/2} (nu - 2)^{1/2} t_nu
nu <- 14
kappa_s_2 <- c(0, 1, 0, 6 / (nu - 4), 0, 240 / ((nu - 4) * (nu - 6)), 0,
               5040 * (5*nu - 22) / ((nu - 4)^2 * (nu - 6) * (nu - 8)))
nu <- 9
kappa_s_3 <- c(0, 1, 0, 6 / (nu - 4), 0, 240 / ((nu - 4) * (nu - 6)), 0,
               5040 * (5*nu - 22) / ((nu - 4)^2 * (nu - 6) * (nu - 8)))


a_s <- c(0, 0.1, 0.3, 0.5, 0.7)
k_s <- 0:2

n_s <- c(25, 50, 75, 100, 150, 200, 250, 500, 1000, 2000)

m_max <- 50

cfg <- list()

for (i_a in 1:length(a_s)) {
  for (i_k in 1:length(k_s)) {
    for (i_dist in 1:3) {
      if (i_dist == 1) {
        cfg <- c(cfg, list(list(a       = a_s[[i_a]],
                                k       = k_s[[i_k]],
                                innov   = c("N(0,1)", "t_14", "t_9")[i_dist],
                                kappa_s = kappa_s_1
        )))
      } else if(i_dist == 2) {
        cfg <- c(cfg, list(list(a       = a_s[[i_a]],
                                k       = k_s[[i_k]],
                                innov   = c("N(0,1)", "t_14", "t_9")[i_dist],
                                kappa_s = kappa_s_2
        )))
      } else if(i_dist == 3) {
        cfg <- c(cfg, list(list(a       = a_s[[i_a]],
                                k       = k_s[[i_k]],
                                innov   = c("N(0,1)", "t_14", "t_9")[i_dist],
                                kappa_s = kappa_s_3
        )))
      }
    }
  }
}

