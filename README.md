# Scope
This repository contains R code that can be used to replicate the results in

> Finite sample distributional error bounds for empirical cross-covariances

There are three parts to this repository:

* Obtaining the Wasserstein distances via simulation,
* computing the values of the bounds, and
* generating figures and tables.

Each part is described separately in the following sections.
The corresponding scripts to the three steps are organised in the folder
`sim` and code that we considered of general use are in the folder `R`.

# Obtaining the Wasserstein distances via simulation
Running the simulations to a high degree of accuracy is computationally
expensive. We have therefore opted for a setup where results can be computed in
parallel. There are 750 tasks in total, obtained by combination of the

* 5 different autoregressive parameter (0, 0.1, 0.3, 0.5, 0.7),
* 3 diffferent distributions for the innovations, and
* B = 50 independent replications.

To replicate the numbers, the following files are relevant:

* `sim/1_Wasserstein_distance/sim.R`,
* `sim/1_Wasserstein_distance/sim_params.R`,
* `R/AR_q.R`.

**Step 1**
The file `sim.R` contains the script to simulate the Wasserstein distances
and save them in files named `W1_[task_id].Rdata`, where `[task_id]` is a
variable that has to be set to take the values between 1-750 every time the
script `sim.R` is run. The 750 different "configurations" are encoded in
the file `sim_params.R`. For running the script on a HPC cluster with SLURM,
an optional shell script named

* `sim/1_Wasserstein_distance/sim.sh`

is provided.

**Step 2**
Next, the simulated Wasserstein distances from step 1, which were stored in
750 separate files are now rearranged into one data.frame `W1`.
To this end, the script `sim/1_Wasserstein_distance/sim_merge.R` has to be run.
The files from step 1 are expected to be in a folder `out/` of the working
directory. This can be changed in line 27 of `sim_merge.R`. The output of step 2,
`W1.Rdata`, is expected to be in `sim/1_Wasserstein_distance` for generation
of the figures and tables (part 3 of the repository).
For convenience the file `W1.Rdata` is also available from the repository.

# Computing the values of the bounds
For the parameters and AR(1) model we consider computing the bound is
computationally less expensive than the simulation task in the previous section.

To replicate the numbers, the following files are relevant:

* `sim/2_bound/sim.R`,
* `sim/2_bound/sim_params.R`,
* `R/indecomposble_partition.R`,
* `R/comp_bound.R`.

Running `sim.R` computes the values of the bound and arranges them in a data
frame `bound_Thm5_AR1` and saves it to the file `bound_Thm5_AR1.Rdata`.
For convenience the file `bound_Thm5_AR1.Rdata` is also available from the
repository.

# Generating figures and tables
The figure and tables in the paper are generated from the output from
part 1 and 2.

To generate the figure and tables call the following scripts:

* `sim/3_tables_figures/figure_1.R`,
* `sim/3_tables_figures/tables.R`.
