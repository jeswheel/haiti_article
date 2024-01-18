# Fitting profiles for Haiti 2.
#
# The fit_haiti2 function is a very simple helper that will take as input
# some starting parameters and an objective function for Model 2 and will
# return parameters estimates along with corresponding estimates of their
# log-likelihood.
#
# This this script we use this helper function to do a profile likelihood search
# of the estimated parameters using the batchtools package. To do this, we will
# create a grid of parameter values and submit jobs where each job will use a
# row from the grid as a starting point to fit the model.


# Load Packages -----------------------------------------------------------

library(batchtools)
library(tidyverse)
library(data.table)
library(haitipkg)
library(pomp)
library(spatPomp)
library(subplex)


# Create Experiment Registry ----------------------------------------------

RUN_LEVEL = 2

reg <- makeExperimentRegistry(
  file.dir = paste0('model2/profileReg_RL', RUN_LEVEL, '_v6'),
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

n_starts <- 1

# Create Profile Design Matrix --------------------------------------------

set.seed(22446688)
h2 <- haiti2(cutoff = 10000, measure = "log")

prof_params <- c(
  "Mu",
  "Beta",
  "BetaW",
  "v",
  "sigma",
  "phase"
)

final_pars <- coef(h2)
prof_vars <- c()
for (pp in prof_params) {

  if (pp == 'Mu') {
    prof_values <- seq(6500, 12500, length.out = 30)
  } else if (pp == "Beta") {
    prof_values <- seq((.Machine$double.eps)^(1/5), (5e-05)^(1/5), length.out = 100)^5
  } else if (pp == "BetaW") {
    prof_values <- seq(0.9, 1.25, length.out = 20)
  } else if (pp == "v") {
    prof_values <- seq(1.2, 1.5, length.out = 20)
  } else if (pp == "sigma") {
    prof_values <- seq((.Machine$double.eps)^(1/5), (1e-03)^(1/5), length.out = 100)^5
  } else if (pp == "phase") {
    prof_values <- seq(6, 8.5, length.out = 30)
  }

  prof_cols <- matrix(rep(prof_values, each = n_starts), ncol = 1)
  colnames(prof_cols) <- pp

  bounds <- tibble::tribble(
    ~param, ~lower, ~upper,
    "Mu",     8000,  10000,
    "Beta",  2e-17,  5e-15,
    "BetaW",  0.75,    1.3,
    "v",      0.75,      2,
    "sigma", 2e-13,  4e-13,
    "phase",  4.36,  10.64
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
    nseq = n_starts * length(prof_values)
  )

  guesses <- dplyr::bind_cols(prof_cols, guesses_tmp)
  all_params <- coef(h2)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]

  fixed_mat <- matrix(
    rep(fixed_params, n_starts * length(prof_values)),
    byrow = TRUE, nrow = n_starts * length(prof_values)
  )

  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(h2))]
  final_pars <- rbind(final_pars, guesses_all)
  prof_vars <- c(prof_vars, rep(pp, nrow(guesses_all)))
}

final_pars <- final_pars[-1, names(coef(h2))]
data_obj <- list(starts = final_pars, prof_vars = prof_vars)


# Define profile problem for registry -------------------------------------

sample_starting_point <- function(
    data = data_obj, job, i, ...
    ) {

  prof_par <- data$prof_vars[i]
  prof_val <- data$starts[i, prof_par]

  all_est_pars <- c('Mu', 'Beta', 'BetaW', 'v', 'sigma', 'phase')
  est_parms <- all_est_pars[all_est_pars != prof_par]
  all_trans_pars <- est_parms[est_parms != 'phase']

  h2 <- haiti2(cutoff = 10000, measure = "log")
  start_parms <- data$starts[i, ]
  start_parms[prof_par] = prof_val

  # Get objective function to minimize
  epi_ofun <- pomp::traj_objfun(
    h2,
    est = est_parms,
    params = start_parms
  )

  # Log-transform all parameters that need to be transformed
  theta <- log(unlist(start_parms[all_trans_pars]))

  # Add "phase" parameter if needed, this variable is not transformed.
  if (prof_par != 'phase') {
    theta <- c(theta, 'phase' = unlist(unname(start_parms['phase'])))
  }

  list(
    'initialization' = theta,
    'obj_fun' = epi_ofun,
    'prof_par' = prof_par,
    'prof_val' = prof_val
  )
}

addProblem(name = 'profile', data = data_obj, fun = sample_starting_point)

# Define algorithm for parameter estimation -------------------------------

fit_model <- function(data, job, instance, ...) {

  # Vector of all parameters that may be fit and transformed
  trans_pars <- c('Mu', 'Beta', 'BetaW', 'v', 'sigma')

  # Get the names of all of the parameters that are fit for this particular
  # instance
  est_pars <- names(instance$initialization)

  # subset the transformed parameter vector using est_pars
  trans_pars <- trans_pars[trans_pars %in% est_pars]

  # Fit the model using the given initialization
  h2_fit <- haitipkg::fit_haiti2(
    initialization = instance$initialization,
    obj_fun = instance$obj_fun
  )

  # Back transform the model parameters to natural scale
  if ('phase' %in% est_pars) {
    estimated_vals <- c(exp(h2_fit$par[trans_pars]), h2_fit$par['phase'])
  } else {
    estimated_vals <- exp(h2_fit$par)
  }

  # Add logLik to result matrix, and rename the final vector
  results <- c('logLik' = -h2_fit$value, estimated_vals, instance$prof_val)
  names(results) <- c('logLik', est_pars, instance$prof_par)

  # Re-order the vector so the results are in the same order for all jobs.
  results[c("logLik", "Mu", "Beta", "BetaW", "v", "sigma", "phase")]
}

addAlgorithm(name = 'fitMod', fun = fit_model)

# Specifying experiment parameters ----------------------------------------

pdes <- list(
  'profile' = data.frame(
    i = 1:nrow(final_pars),
    prof_var = prof_vars
  )
)

addExperiments(prob.designs = pdes)

# Submit Jobs -------------------------------------------------------------

resources1 <- list(
  account = 'stats_dept1', walltime = '45:00',
  memory = '1000m', ncpus = 1
)

submitJobs(
  data.table(
    job.id = 1:nrow(final_pars),
    chunk = 1:(round(nrow(final_pars)) / 3)
  ), resources = resources1
)

