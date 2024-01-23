# Fitting profiles for Haiti 1.
#
# ********
# NOTE:
#
#   If you would like to run this code but have never used the `batchtools` R
#   package, then we strongly recommend completing both Examples 1 and 2 from
#   the batchtools vignette, by first installing the `batchtools` package and
#   running: `vignette('batchtools')`, or visiting:
#
#         https://mllg.github.io/batchtools/articles/batchtools.html
#
# WARNING:
#
#   This script was ran using a high-performance computing cluster. In total,
#   the RUN_LEVEL == 3 version of this script took TODO compute hours. Partial
#   results can be obtained by lowering the required computational effort
#   (params: NP, NMIF, nprof, NP_EVAL, NREPS_EVAL), which can be done by
#   reducing the RUN_LEVEL to 2 or 1, or alternatively by manually changing
#   these values.
# ********
#
# The haitipkg::fit_haiti1 function is a very simple helper that will take as
# input some starting parameters, and hyperparameters for the IF2 algorithm, and
# will return parameter estimates along with corresponding estimates of their
# log-likelihoods.
#
# In this script we will use this function to do a profile likelihood search
# of all model parameters using the `batchtools` package. To do this, we will
# create a grid of parameter values, and submit jobs where each job will use
# a row from the grid as a starting point to fit the model.
#
# The number of rows in the grid (which is the `final_pars` object that is
# created using a for-loop over the profile parameters) correspond to the total
# number of jobs that will be submitted using `batchtools`. Each job takes
# roughly the same amount of time to compute, and hence the run-time complexity
# is linear with respect to the total number of jobs. The total number of jobs
# can be controlled by either (1) changing the number of variables to profile
# over (2) changing the number of unique values for each profile parameter (3)
# changing the number of points per profile point.
#
# The time-complexity of each job is roughly linear in the parameters NP, and
# NMIF. The fitting algorithm (IF2) has a computational complexity that is
# linear with these parameters, and the fitting procedure takes the majority
# of the time for each job. The model fit is evaluated in each job using the
# particle filter algorithm, which has a computational complexity that is
# linear in NP_EVAL and NREPS_EVAL, so modifying these parameters will also
# have an effect on the total computational effort required, but is typically
# much less than the effort used to fit the model.

# Load Packages -----------------------------------------------------------

library(batchtools)
library(tidyverse)
library(data.table)
library(haitipkg)

RUN_LEVEL <- 3

nprof      <- switch(RUN_LEVEL,  2,    15, 25)
NMIF       <- switch(RUN_LEVEL,  5,   150, 200)
NP         <- switch(RUN_LEVEL,  50, 1000, 1500)
NP_EVAL    <- switch(RUN_LEVEL, 100, 3000, 5000)
NREPS_EVAL <- switch(RUN_LEVEL,   3,   10, 15)
COOLING    <- 0.5

# Create Experiment Registry ----------------------------------------------

reg <- makeExperimentRegistry(
  file.dir = paste0('model1/profileReg_RL', RUN_LEVEL),
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

# Create Profile Design Matrix --------------------------------------------

set.seed(665544)

# Create Model 1 pomp object, used to get fixed parameters in the searches.
h1 <- haiti1_joint()

# Vector of parameters to profile over.
prof_params <- c(
  'betat',
  'tau_epi',
  'tau_end',
  'rho',
  'nu',
  'sig_sq_epi',
  'sig_sq_end',
  'E_0',
  'I_0',
  'beta1',
  'beta2',
  'beta3',
  'beta4',
  'beta5',
  'beta6'
)

# Initialize parameter matrix using default parameter values.
final_pars <- coef(h1)
prof_vars <- c()
for (pp in prof_params) {  # Loop through all profile parameters

  # Profile bounds estimated using RL = 2
  if (pp == "betat") {
    prof_values <- seq(-0.15, 0.05, length.out = 30)
  } else if (pp == 'tau_epi') {
    prof_values <- c(exp(seq(2.75, 9, length.out = 40)), 5e+12)
  } else if (pp == 'tau_end') {
    prof_values <- c(exp(seq(2, 9, length.out = 40)), 5e+12)
  } else if (pp == 'rho') {
    prof_values <- seq(0.1, 0.9, length.out = 30)
  } else if (pp == 'nu') {
    prof_values <- seq(0.9, 1, length.out = 15)
  } else if (pp == 'sig_sq_epi') {
    prof_values <- seq(0.04, 0.2, length.out = 30)
  } else if (pp == 'sig_sq_end') {
    prof_values <- seq(0.02, 0.275, length.out = 30)
  } else if (pp == 'E_0') {
    prof_values <- seq(1 / 10911819, 25000 / 10911819, 1000 / 10911819)
  } else if (pp == 'I_0') {
    prof_values <- seq(1 / 10911819, 25000 / 10911819, 1000 / 10911819)
  } else if (pp == 'beta1') {
    prof_values <- seq(0.5, 2.5, length.out = 20)
  } else if (pp == 'beta2') {
    prof_values <- seq(0.25, 2, length.out = 20)
  } else if (pp == 'beta3') {
    prof_values <- seq(0.4, 2.1, length.out = 20)
  } else if (pp == 'beta4') {
    prof_values <- seq(0.4, 2.1, length.out = 20)
  } else if (pp == 'beta5') {
    prof_values <- seq(0.5, 2.5, length.out = 20)
  } else if (pp == 'beta6') {
    prof_values <- seq(-0.1, 2, length.out = 20)
  }

  prof_cols <- matrix(rep(prof_values, each = nprof), ncol = 1)
  colnames(prof_cols) <- pp

  # Parameter bounds found using RL 2.
  bounds <- tibble::tribble(
    ~param, ~lower, ~upper,
    "beta1",        1.00,  2.15,
    "beta2",         .50,  2.00,
    "beta3",         .75,  1.75,
    "beta4",         .75,  1.50,
    "beta5",        1.00,  2.10,
    "beta6",         .50,  1.50,
    "betat",       -0.15,  0.05,
    "tau_epi",       180,  1250,
    "tau_end",        50,  1000,
    "rho",          0.2,   0.95,
    "nu",           0.94,     1,
    "sig_sq_epi",  0.075, 0.125,
    "sig_sq_end",   0.05,   0.2,
    "E_0", 10 / 10911819, 5000 / 10911819,
    "I_0", 10 / 10911819, 20000 / 10911819
  )

  bounds <- bounds %>%
    filter(param != pp)

  lower <- bounds$lower
  names(lower) <- bounds$param

  upper <- bounds$upper
  names(upper) <- bounds$param

  guesses_tmp <- pomp::runif_design(
    lower = lower,
    upper = upper,
    nseq = nprof * length(prof_values)
  )

  guesses <- dplyr::bind_cols(prof_cols, guesses_tmp)

  S_0 <- matrix(unname(1 - guesses[, 'I_0'] - guesses[, 'E_0']), ncol = 1)
  colnames(S_0) <- 'S_0'
  guesses <- dplyr::bind_cols(guesses, S_0)

  # We need to add fixed parameters
  all_params <- coef(h1)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]

  fixed_mat <- matrix(
    rep(fixed_params, nprof * length(prof_values)),
    byrow = TRUE, nrow = nprof * length(prof_values)
  )

  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(h1))]
  final_pars <- rbind(final_pars, guesses_all)
  prof_vars <- c(prof_vars, rep(pp, nrow(guesses_all)))
}

final_pars <- final_pars[-1, names(coef(h1))]

data_obj <- list(starts = final_pars, prof_vars = prof_vars)
# dim(final_pars)

# Define profile problem for registry -------------------------------------

# A simple function that will sample the starting values of the profile search
sample_starting_point <- function(data = data_obj, job, i, ...) {
  list(
    start_par = unlist(data$starts[i, ]),
    prof_var = data$prof_vars[i]
  )
}

addProblem(name = 'profile', data = data_obj, fun = sample_starting_point)

# Define algorithm for parameter estimation -------------------------------

# This function estimates model parameters using the IF2 algorithm. All
# calibrated parameters besides the profile parameter in the instance object
# and the fixed parameters have a non-zero RW_SD.
fit_model <- function(
    data, job, instance, nmif, np, np_eval, nreps_eval,
    cooling, ...
) {

  est_param_names <- c(
    paste0('beta', 1:6), 'betat', 'rho', 'nu', 'E_0', 'I_0',
    'tau_epi', 'tau_end', 'sig_sq_epi', 'sig_sq_end'
  )

  est_param_names <- est_param_names[est_param_names != instance$prof_var]

  # Simple function that is used to set the rw_sd.
  # This function is convenient so that we don't have to type out the rw_sd
  # for each unit parameter.
  set_rw <- function(x) {
    if (x == 'nu') {
      return(0.01)
    } else if (x == 'E_0' || x == "I_0") {
      return(expression(ivp(0.2)))
    } else if (x == 'tau_epi' || x == 'sig_sq_epi') {
      return(expression(ifelse(time < 232, 0.02, 0)))
    } else if (x == 'tau_end' || x == 'sig_sq_end') {
      return(expression(ifelse(time >= 232, 0.02, 0)))
    } else {
      return(0.02)
    }
  }

  reg_rw.sd <- lapply(est_param_names, set_rw)
  names(reg_rw.sd) <- est_param_names

  RW_SD <- do.call(pomp::rw_sd, reg_rw.sd)

  fit_haiti1(
    start_params = instance$start_par,
    NMIF = nmif,
    NP = np,
    NP_EVAL = np_eval,
    NREPS_EVAL = nreps_eval,
    RW_SD = RW_SD,
    COOLING = cooling
  )
}

addAlgorithm(name = 'fitMod', fun = fit_model)

# Completing the experiment -----------------------------------------------

pdes <- list('profile' = data.frame(i = 1:nrow(final_pars), prof_var = prof_vars))
ades <- list(
  'fitMod' = data.frame(
    nmif = NMIF, np = NP, np_eval = NP_EVAL, nreps_eval = NREPS_EVAL, cooling = COOLING
  )
)

addExperiments(prob.designs = pdes, algo.designs = ades)

# Submit Jobs -------------------------------------------------------------

# resources1 <- list(account = 'stats_dept1', walltime = '10:00', memory = '1000m', ncpus = 1)
resources1 <- list(
  account = 'stats_dept1', walltime = '3:30:00',
  memory = '5000m', ncpus = 1
  )

submitJobs(
   data.table(
      job.id = 1:nrow(final_pars),
      chunk = 1:(round(nrow(final_pars) / 3))
   ), resources = resources1
)

# Get Results -------------------------------------------------------------

# You can check the status of the jobs (how many are queued, how many have
# started, how many have finished, and how many have errors / ran out of time)
# by running:
#
# getJobStatus()
#
# You can ensure that you won't summarize the results until everything has
# finished by running:
#
# waitForJobs()
#
# Because there are many long jobs, however, it may not be the best idea to do
# this interactively. Instead, it's wise to periodically check on the status
# of the jobs when you think they may have finished. Therefore if you would
# like to check on your results after closing your R session, you need to
# re-load your registry (after opening up a new R-session in the same
# working directory that was used to create the registry) using:
#
# library(batchtools)
#
# RUN_LEVEL = 3  # Or adjust RUN_LEVEL as needed.
# reg = loadRegistry(paste0('model1/profileReg_RL', RUN_LEVEL))
# getStatus()
#
# Once all of the jobs are finished, we can summarize the results using:
#
# h1_profile_results = unwrap(reduceResultsDataTable())
# h1_pars = unwrap(getJobPars())
# h1_results = ijoin(h1_pars, h1_profile_results)
#
# saveRDS(
#   h1_profile_results,
#   paste0('model1/run_level_', RUN_LEVEL, '/', 'h1_profiles.rds')
# )
#
# We can also get information about each job, which can be helpful to pick the
# size of each run level and to estimate the cost of the RUN_LEVEL_3 search.
# For example, we can see the total / average computational time of all of the
# jobs by running:
#
# h1_stats <- unwrap(getJobStatus())
# h1_stats <- ijoin(h1_pars, h1_stats)
#
# summary(h1_stats$time.running |> as.numeric())  # unit of measurement is seconds
