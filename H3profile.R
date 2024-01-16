# Fitting profiles for Haiti 3.
#
# The fit_haiti3 function is a very simple helper that will take as input
# some starting parameters, and hyperparameters for the IBPF algorithm, and
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

nprof <- switch(RUN_LEVEL, 2, 14, 20)
NBPF <- switch(RUN_LEVEL, 5, 50, 100)
NP <- switch(RUN_LEVEL, 50, 750, 1000)
NP_EVAL <- switch(RUN_LEVEL, 100, 1000, 2000)
NREPS_EVAL <- switch(RUN_LEVEL, 3, 6, 10)
SPAT_REGRESSION <- 0.05
COOLING <- 0.3

# Create Experiment Registry ----------------------------------------------

reg <- makeExperimentRegistry(
  file.dir = paste0('model3/profileReg_RL', RUN_LEVEL, '_v6'),
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

# Create Profile Design Matrix --------------------------------------------

set.seed(665544)

# Read in previous best results
haiti3_fit <- readRDS("model3/run_level_3/haiti3_fit.rds")
best_pars <- haiti3_fit$search2$params[which.max(haiti3_fit$search2$logLiks$logLik), ]
h3_spat <- haiti3_spatPomp()

# Create vectors for the unit and shared parameters
unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")
shared_param_names <- c(
  "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
  "epsilon", "k"
)

est_param_names <- c(
  unit_specific_names, shared_param_names
)

# prof_params <- shared_param_names
prof_params <- c(
  "XthetaA"
)

final_pars <- best_pars
prof_vars <- c()
for (pp in prof_params) {

  if (pp == "mu_B") {
    prof_values <- seq(200, 650, length.out = 25)
  } else if (pp == 'XthetaA') {
    prof_values <- seq(0.002, 0.125, length.out = 31)
  } else if (pp == 'thetaI') {
    prof_values <- seq(0.000075, 0.00016, length.out = 20)
  } else if (pp == 'lambdaR') {
    prof_values <- seq(1, 3, length.out = 20)
  } else if (pp == 'r') {
    prof_values <- seq(0.15, 1, length.out = 21)
  } else if (pp == 'std_W') {
    prof_values <- seq(0.025, 0.05, length.out = 25)
  } else if (pp == 'epsilon') {
    prof_values <- seq(0.85, 0.9999, length.out = 20)
  } else if (pp == 'k') {
    prof_values <- seq(10, 95, length.out = 25)
  }

  tmp_pars <- matrix(
    rep(best_pars, length(prof_values)),
    nrow = length(prof_values),
    ncol = length(best_pars),
    byrow = TRUE
  )

  colnames(tmp_pars) <- names(best_pars)

  if (pp %in% shared_param_names) {
    prof_cols <- matrix(rep(rep(prof_values, 10), each = nprof), ncol = 10)
    colnames(prof_cols) <- paste0(pp, 1:10)
    tmp_pars[, paste0(pp, 1:10)] <- prof_values
  } else {
    prof_cols <- matrix(rep(prof_values, each = nprof), ncol = 1)
    colnames(prof_cols) <- pp
    tmp_pars[, pp] <- prof_values
  }

  # Unit-lower bounds
  unit_lower_bounds <- c(
    'betaB1' = 2.6,
    'betaB2' = 9,
    'betaB3' = 13,
    'betaB4' = 13,
    'betaB5' = 2.5,
    'betaB6' = 15,
    'betaB7' = 5,
    'betaB8' = 0.4,
    'betaB9' = 5,
    'betaB10' = 6,
    'foi_add1' = 5e-8,
    'foi_add2' = 1e-8,
    'foi_add3' = 1e-7,
    'foi_add4' = 5e-8,
    'foi_add5' = 2e-7,
    'foi_add6' = 2.5e-7,
    'foi_add7' = 1.1e-7,
    'foi_add8' = 1e-8,
    'foi_add9' = 1.5e-7,
    'foi_add10' = 1e-7,
    'aHur3' = 0.01,
    'aHur9' = 10,
    'hHur3' = 30,
    'hHur9' = 40,
    'Iinit3' = 0.9 / 468301,
    'Iinit4' = 4.8 / 342525
  )

  # Unit upper-bounds
  unit_upper_bounds <- c(
    'betaB1' = 8,
    'betaB2' = 35,
    'betaB3' = 50,
    'betaB4' = 50,
    'betaB5' = 10,
    'betaB6' = 50,
    'betaB7' = 20,
    'betaB8' = 2,
    'betaB9' = 20,
    'betaB10' = 20,
    'foi_add1' = 8e-7,
    'foi_add2' = 6e-7,
    'foi_add3' = 5e-7,
    'foi_add4' = 2.5e-7,
    'foi_add5' = 7.5e-7,
    'foi_add6' = 6e-7,
    'foi_add7' = 4.5e-7,
    'foi_add8' = 4e-7,
    'foi_add9' = 3.5e-7,
    'foi_add10' = 2.8e-7,
    'aHur3' = 50,
    'aHur9' = 45,
    'hHur3' = 120,
    'hHur9' = 100,
    'Iinit3' = 35 / 468301,
    'Iinit4' = 35 / 342525
  )

  shared_lower_bounds <- c(
    "mu_B" = 450,
    "XthetaA" = 0.01,
    "thetaI" = 2.5e-05,
    "lambdaR" = 0.75,
    "r" = 0.2,
    "std_W" = 0.0275,
    "epsilon" = 0.8,
    "k" = 25
  )

  shared_upper_bounds <- c(
    "mu_B" = 700,
    "XthetaA" = 0.15,
    "thetaI" = 1e-04,
    "lambdaR" = 3.5,
    "r" = 1,
    "std_W" = 0.0375,
    "epsilon" = 0.99,
    "k" = 80
  )

  est_u_names <- names(unit_lower_bounds)
  partial_u_names <- est_u_names[est_u_names != pp]
  partial_s_names <- shared_param_names[shared_param_names != pp]

  unit_lower_bounds <- unit_lower_bounds[partial_u_names]
  unit_upper_bounds <- unit_upper_bounds[partial_u_names]

  shared_lower_bounds <- shared_lower_bounds[partial_s_names]
  shared_upper_bounds <- shared_upper_bounds[partial_s_names]

  # Create data.frame with random unit parameters
  guesses_unit <- pomp::runif_design(
    lower = unit_lower_bounds,
    upper = unit_upper_bounds,
    nseq = nprof * length(prof_values)
  )

  guesses_shared <- pomp::runif_design(
    lower = shared_lower_bounds,
    upper = shared_upper_bounds,
    nseq = nprof * length(prof_values)
  )

  # Need to duplicate each of the shared parameter columns
  guesses_shared <- guesses_shared[, rep(partial_s_names, each = 10)]
  colnames(guesses_shared) <- paste0(rep(partial_s_names, each = 10), 1:10)

  # Combine the unit and shared parameters
  guesses <- cbind(guesses_shared, guesses_unit)
  guesses <- dplyr::bind_cols(prof_cols, guesses)

  # We need to add fixed parameters
  all_params <- coef(h3_spat)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  fixed_mat <- matrix(
    rep(fixed_params, nprof * length(prof_values)),
    byrow = TRUE, nrow = nprof * length(prof_values)
  )
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(h3_spat))]
  guesses_all <- rbind(guesses_all, tmp_pars)[, names(coef(h3_spat))]
  final_pars <- rbind(final_pars, guesses_all)
  prof_vars <- c(prof_vars, rep(pp, nrow(guesses_all)))
}

final_pars <- final_pars[-1, names(coef(h3_spat))]

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
    data, job, instance, nbpf, np, spat_regression, np_eval, nreps_eval,
    cooling, start_date = "2010-11-20", ...
) {

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")
  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "std_W",
    "epsilon", "k"
  )

  est_param_names <- c(
    unit_specific_names, shared_param_names
  )

  # Add unit numbers to each parameter
  est_param_names_expanded <- paste0(rep(est_param_names, each = 10), 1:10)

  # Simple function that is used to set the rw_sd for the first search
  # This function is convenient so that we don't have to type out the rw_sd
  # for each unit parameter.
  set_rw_1 <- function(x) {

    if (x == instance$prof_var || gsub("[[:digit:]]+$", "", x) == instance$prof_var) {
      return(0)
    } else if (gsub("[[:digit:]]+$", "", x) %in% shared_param_names) {
      return(0.005)
    } else if (x %in% paste0(rep(c("aHur", 'hHur'), each = 2), c(3, 9))) {
      return(expression(ifelse(time >= 2016.754 & time <= 2017, 0.01, 0)))
    } else if (x %in% c("Iinit3", "Iinit4")) {
      return(expression(ivp(0.15)))
    } else if (grepl('^foi_add[[:digit:]]+$', x)) {
      return(0.0075)
    } else if (grepl('^betaB[[:digit:]]+$', x)) {
      return(0.005)
    } else {
      return(0)
    }
  }

  reg_rw_1.sd <- lapply(est_param_names_expanded, set_rw_1)
  names(reg_rw_1.sd) <- est_param_names_expanded

  RW_SD <- do.call(rw_sd, reg_rw_1.sd)

  fit_haiti3(
    start_params = instance$start_par,
    NBPF = nbpf,
    NP = np,
    SPAT_REGRESSION = spat_regression,
    NP_EVAL = np_eval,
    NREPS_EVAL = nreps_eval,
    RW_SD = RW_SD,
    COOLING = cooling,
    start_date = start_date
  )
}

addAlgorithm(name = 'fitMod', fun = fit_model)

# Completing the experiment -----------------------------------------------

pdes <- list('profile' = data.frame(i = 1:nrow(final_pars)))
ades <- list(
  'fitMod' = data.frame(
    nbpf = NBPF, np = NP, spat_regression = SPAT_REGRESSION,
    np_eval = NP_EVAL, nreps_eval = NREPS_EVAL, cooling = COOLING
  )
)

addExperiments(prob.designs = pdes, algo.designs = ades)

# Submit Jobs -------------------------------------------------------------

# resources1 <- list(account = 'stats_dept1', walltime = '10:00', memory = '5000m', ncpus = 1)
resources1 <- list(account = 'stats_dept1', walltime = '2:30:00', memory = '5000m', ncpus = 1)

submitJobs(resources = resources1)
# submitJobs(ids = 1, resources = resources1)
