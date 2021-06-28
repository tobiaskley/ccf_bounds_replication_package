# This script computes the bound for the AR(1) model configurations in
# sim_params.R.
# 
# Author: KLEY
###############################################################################

library(parallel)

source("sim_params.R")
source("../../R/indecomposablePartitions.R")
source("../../R/comp_bound.R")

# data frames for results
out <- data.frame(a = numeric(),
                  n = numeric(),
                  k = numeric(),
                  innov = character(),
                  bound_star = numeric(),
                  m_star = numeric())

out$innov <- factor(out$innov, levels = c(1, 2, 3),
                    labels = c("N(0,1)", "t_14", "t_9"))

for (task_id in 1:length(cfg)) {

  # detach parameters from configuration
  a       <- cfg[[task_id]]$a
  k       <- cfg[[task_id]]$k
  innov   <- cfg[[task_id]]$innov
  kappa_s <- cfg[[task_id]]$kappa_s

  res <- comp_bound_AR1(n_s, a, k, kappa_s, 1, m_max)
  for (i_n in 1:length(n_s)) {
    
    out <- rbind(out, data.frame(a = cfg[[task_id]]$a,
                                 n = n_s[i_n],
                                 k = cfg[[task_id]]$k,
                                 innov = cfg[[task_id]]$innov,
                                 bound_star = res$results[[i_n]]$bound_star,
                                 m_star = res$results[[i_n]]$m_star,
                                 m_max = max(res$m_s),
                                 is_local_min = (min(res$results[[i_n]]$bound_s) <
                                                   tail(res$results[[i_n]]$bound_s, 1))
    ))
  }

}

bound_Thm5_AR1 <- out
save(bound_Thm5_AR1, file = "bound_Thm5_AR1.Rdata")
