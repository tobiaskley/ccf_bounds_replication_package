# This script merges the results computed on the cluster to one data frame
# and saves it as W1.Rdata
# 
# Author: KLEY
###############################################################################

library(parallel)

source("../../R/AR_q.R")
source("sim_params.R")

# data frames for results
W1 <- data.frame(cfg_id = numeric(),
                  a = numeric(),
                  innov = character(),
                  b = numeric(),
                  n = numeric(),
                  k = numeric(),
                  W1 = numeric())

W1$innov <- factor(W1$innov, levels = 1:length(innov_name_s), labels = innov_name_s)

files_completed <- 0
for (task_id in 1:length(cfg)) {

  # check if results are available
  file_name <- paste("out/W1_", task_id, ".Rdata", sep = "")

  if (file.exists(file_name)) {
    load(file_name)

    for (b in b_start:b_end) {
      for (i_n in 1:length(n_s)) {
        for (i_k in 1:length(k_s)) {
          W1 <- rbind(W1, data.frame(cfg_id = task_id,
                                     a = a,
                                     innov = innov_name,
                                     b = b,
                                     n = n_s[i_n],
                                     k = k_s[i_k],
                                     W1 = W1_gamma[b - b_start + 1, i_n, i_k]
                                     ))
        }
      }
    }

    cat("task id: ", task_id, ", loaded ", file_name, "\n")
    files_completed <- files_completed + 1

  }
}

save(W1, file = "W1.Rdata")
