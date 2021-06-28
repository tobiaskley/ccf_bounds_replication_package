# This script returns Tables 1-6.
# 
# Author: KLEY
###############################################################################


load("../1_Wasserstein_distance/W1.Rdata")
load("../2_bound/bound_Thm5_AR1.Rdata")

for (i in 1:nrow(bound_Thm5_AR1)) {
  
  bound_Thm5_AR1$W1[i] <- mean(W1[W1$innov == bound_Thm5_AR1$innov[i] &
                                       W1$n == bound_Thm5_AR1$n[i] &
                                       W1$k == bound_Thm5_AR1$k[i] &
                                       W1$a == bound_Thm5_AR1$a[i], ]$W1)
  
  bound_Thm5_AR1$W1_sd[i] <- sd(W1[W1$innov == bound_Thm5_AR1$innov[i] &
                                     W1$n == bound_Thm5_AR1$n[i] &
                                     W1$k == bound_Thm5_AR1$k[i] &
                                     W1$a == bound_Thm5_AR1$a[i], ]$W1)
  
}


bound_Thm5_AR1$coef_var <- bound_Thm5_AR1$W1_sd / bound_Thm5_AR1$W1
bound_Thm5_AR1$ratio <- bound_Thm5_AR1$bound_star / bound_Thm5_AR1$W1

# next two lines can be adapted to get smaller tables
bound_Thm5_AR1 <- bound_Thm5_AR1[bound_Thm5_AR1$a %in% c(0, 0.1, 0.3, 0.5, 0.7),]
bound_Thm5_AR1 <- bound_Thm5_AR1[bound_Thm5_AR1$k %in% 0:2,]

write.csv(bound_Thm5_AR1, file = "bound_AR1_example.csv")

library(xtable)

# bounds as a table
for (i_innov in 1:3) {
  
  innov <- c("N(0,1)", "t_14", "t_9")[i_innov]
  innov_latex <- c("\\mathcal{N}(0,1)", "\\sqrt{12/14} \\ t_{14}", "\\sqrt{7/9} \\ t_{9}")[i_innov]
  innov_label <- c("N01", "t14", "t9")[i_innov]
  
  xt <- xtabs(round(bound_star, 3) ~ a + n + k,
              data = bound_Thm5_AR1[bound_Thm5_AR1$a >= 0 &
                                      bound_Thm5_AR1$innov == innov,])
  ft <- ftable(xt, row.vars = c("k", "a"), col.vars = c("n"))
  
  print(xtableFtable(ft, method = "compact",
                     caption = paste("Value of the bound from Theorem~\\ref{general_theorem} in combination with \\eqref{bnd:Smn_new_bound}, ",
                               "with $m = m^*$ to minimise the bound as described in Section~\\ref{sec:comp_bnd}, ",
                               "for empirical autocovariances, for a range of lags $k$ and sample sizes $n$. ",
                               "The data stems from an AR(1) process with $\\varepsilon_t \\sim ", innov_latex, "$ where $a$ takes a range of values. ",
                               "\\label{tab:bnd_", innov_label,"}", sep=""),
                     table.placement = "ht",
                     size = "\\footnotesize"),
        file = paste("bound_",innov ,".tex", sep = "")
  )
  
  xt <- xtabs(round(W1, 3) ~ a + n + k,
              data = bound_Thm5_AR1[bound_Thm5_AR1$a >= 0 &
                                      bound_Thm5_AR1$innov == innov,])
  ft <- ftable(xt, row.vars = c("k", "a"), col.vars = c("n"))
  
  print(xtableFtable(ft, method = "compact",
                     caption = paste("Value of the true 1-Wasserstein distance considered in Theorem~\\ref{general_theorem} for empirical autocovariances, for ",
                                     "a range of lags $k$ and sample sizes $n$. ",
                                     "The data stems from an AR(1) process with $\\varepsilon_t \\sim ", innov_latex, "$ where $a$ takes a range of values. ",
                                     "\\label{tab:W1_", innov_label,"}", sep=""),
                     table.placement = "ht",
                     size = "\\footnotesize"),
        file = paste("W1_",innov ,".tex", sep = "")
  )
  
}

