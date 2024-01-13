# Fitting profiles for Haiti 1.
#
# The fit_haiti1 function is a very simple helper that will take as input
# some starting parameters, and hyperparameters for the IF2 algorithm, and
# will return parameter estimates along with corresponding estimates of their
# log-likelihoods.
#
# In this script we will use this function to do a profile likelihood search
# of all model parameters using the batchtools package. To do this, we will
# create a grid of parameter values, and submit jobs where each job will use
# a row from the grid as a starting point to fit the model.


# Load Packages -----------------------------------------------------------

library(batchtools)
library(tidyverse)
# library(data.table)
library(haitipkg)

RUN_LEVEL <- 2

nprof      <- switch(RUN_LEVEL,  2,    14, 20)
NMIF       <- switch(RUN_LEVEL,  5,   100, 200)
NP         <- switch(RUN_LEVEL,  50, 1000, 1000)
NP_EVAL    <- switch(RUN_LEVEL, 100, 2500, 10000)
NREPS_EVAL <- switch(RUN_LEVEL,   3,   10, 10)
COOLING    <- 0.5

# Create Experiment Registry ----------------------------------------------

reg <- makeExperimentRegistry(
  file.dir = paste0('model1/profileReg_RL', RUN_LEVEL, '_v1'),
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

# Create Profile Design Matrix --------------------------------------------

set.seed(665544)

# Read in previous best results
haiti1_fit <- readRDS("model1/run_level_3/haiti1_fit.rds")
best_pars <- haiti1_fit %>%
  filter(ll == max(ll)) %>%
  select(-ll, -ll.se) %>%
  unlist()

h1 <- haiti1_joint()
coef(h1) <- best_pars

prof_params <- c(
  'betat',
  'tau_epi',
  'tau_end',
  'rho',
  'nu',
  'sig_sq_epi',
  'sig_sq_end',
  'E_0',
  'I_0'
)

final_pars <- best_pars
prof_vars <- c()
for (pp in prof_params) {

  if (pp == "betat") {
    prof_values <- seq(-0.15, 0.05, length.out = 30)
  } else if (pp == 'tau_epi') {
    prof_values <- seq(180, 1800, length.out = 30)
  } else if (pp == 'tau_end') {
    prof_values <- seq(50, 1800, length.out = 30)
  } else if (pp == 'rho') {
    prof_values <- seq(0.1, 1, length.out = 30)
  } else if (pp == 'nu') {
    prof_values <- seq(0.94, 1, length.out = 20)
  } else if (pp == 'sig_sq_epi') {
    prof_values <- seq(0.075, 0.125, length.out = 30)
  } else if (pp == 'sig_sq_end') {
    prof_values <- seq(0.05, 0.2, length.out = 30)
  } else if (pp == 'E_0') {
    prof_values <- seq(10 / 10911819, 5000 / 10911819, 100 / 10911819)
  } else if (pp == 'I_0') {
    prof_values <- seq(10 / 10911819, 25000 / 10911819, 500 / 10911819)
  }

  tmp_pars <- matrix(
    rep(best_pars, length(prof_values)),
    nrow = length(prof_values),
    ncol = length(best_pars),
    byrow = TRUE
  )

  colnames(tmp_pars) <- names(best_pars)

  prof_cols <- matrix(rep(prof_values, each = nprof), ncol = 1)
  colnames(prof_cols) <- pp
  tmp_pars[, pp] <- prof_values

  bounds <- tibble::tribble(
    ~param, ~lower, ~upper,
    "beta1",        1.00,  2.15,
    "beta2",         .75,  1.75,
    "beta3",         .75,  1.75,
    "beta4",         .75,  1.75,
    "beta5",        1.00,  2.10,
    "beta6",         .50,  1.50,
    "betat",       -0.15,  0.05,
    "tau_epi",       180,  1800,
    "tau_end",        50,  1800,
    "rho",          0.15,     1,
    "nu",           0.95,     1,
    "sig_sq_epi",  0.075, 0.125,
    "sig_sq_end",   0.05,   0.2,
    "E_0", 10 / 10911819, 5000 / 10911819,
    "I_0", 10 / 10911819, 25000 / 10911819
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
  guesses_all <- rbind(guesses_all, tmp_pars)[, names(coef(h1))]
  final_pars <- rbind(final_pars, guesses_all)
  prof_vars <- c(prof_vars, rep(pp, nrow(guesses_all)))
}

final_pars <- final_pars[-1, names(coef(h1))]

data_obj <- list(starts = final_pars, prof_vars = prof_vars)
# dim(final_pars)

# Define profile problem for registry -------------------------------------

sample_starting_point <- function(data = data_obj, job, i, ...) {
  list(
    start_par = unlist(data$starts[i, ]),
    prof_var = data$prof_vars[i]
  )
}

addProblem(name = 'profile', data = data_obj, fun = sample_starting_point)

# Define algorithm for parameter estimation -------------------------------

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

pdes <- list('profile' = data.frame(i = 1:nrow(final_pars)))
ades <- list(
  'fitMod' = data.frame(
    nmif = NMIF, np = NP, np_eval = NP_EVAL, nreps_eval = NREPS_EVAL, cooling = COOLING
  )
)

addExperiments(prob.designs = pdes, algo.designs = ades)

# Submit Jobs -------------------------------------------------------------

# resources1 <- list(account = 'stats_dept1', walltime = '10:00', memory = '1000m', ncpus = 1)
resources1 <- list(account = 'stats_dept1', walltime = '50:00', memory = '5000m', ncpus = 1)

submitJobs(
   data.table(
      job.id = 1:nrow(final_pars), 
      chunk = 1:(round(nrow(final_pars) / 2))
   ), resources = resources1
)
# submitJobs(ids = 1, resources = resources1)
