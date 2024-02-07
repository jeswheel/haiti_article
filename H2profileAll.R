# Fitting profiles for Haiti 2.
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
# ********
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

reg <- makeExperimentRegistry(
  file.dir ='model2/profileAll',
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

# More than one start could be considered, but our preliminary results 
# suggest that only one is needed. 
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
  "phase",
  paste0("InitInfected", 1:10),
  'Delta',
  'Rho', 
  'redmu', 
  'AlphaS'
)

final_pars <- coef(h2)
prof_vars <- c()
for (pp in prof_params) {

  if (pp == 'Mu') {
    prof_values <- seq(6750, 13000, length.out = 30)
  } else if (pp == "Beta") {
    prof_values <- seq((.Machine$double.eps)^(1/5), (6.5e-05)^(1/5), length.out = 50)^5
  } else if (pp == "BetaW") {
    prof_values <- seq(1, 1.2, length.out = 30)
  } else if (pp == "v") {
    prof_values <- seq(1.225, 1.44, length.out = 25)
  } else if (pp == "sigma") {
    prof_values <- seq((.Machine$double.eps)^(1/5), (2e-03)^(1/5), length.out = 50)^5
  } else if (pp == "phase") {
    prof_values <- seq(6.5, 7.9, length.out = 30)
  } else if (pp == 'InitInfected1') {
    prof_values <- seq(5000, 40000, length.out = 50)
  } else if (pp == 'InitInfected2') {
    prof_values <- seq(1000, 7000, length.out = 40)
  } else if (pp == "initInfected 3") {
    prof_values <- seq(0.5, 1000, length.out = 30)
  } else if (pp == "initInfected 4") {
    prof_values <- seq(0.5, 1000, length.out = 30)
  } else if (pp == "initInfected 5") {
    prof_values <- seq(0.5, 1500, length.out = 30)
  } else if (pp == "initInfected 6") {
    prof_values <- seq(0.5, 1000, length.out = 30)
  } else if (pp == "initInfected 7") {
    prof_values <- seq(0.5, 1000, length.out = 30)
  } else if (pp == "initInfected 8") {
    prof_values <- seq(1000, 5000, length.out = 40)
  } else if (pp == "initInfected 9") {
    prof_values <- seq(0.5, 1000, length.out = 30)
  } else if (pp == "initInfected 10") {
    prof_values <- seq(0.5, 1000, length.out = 30)
  } else if (pp == 'Delta') {
    prof_values <- seq(1, 30, 1)
  } else if (pp == 'rho') {
    prof_values <- seq(0.4, 0.99, length.out = 30)
  } else if (pp == 'redmu') {
    prof_values <- seq(5e-8, 5e-5, length.out = 35)
  } else if (pp == 'AlphaS') {
    prof_values <- seq(0.2, 0.6, length.out = 20)
  }

  prof_cols <- matrix(rep(prof_values, each = n_starts), ncol = 1)
  colnames(prof_cols) <- pp

  bounds <- tibble::tribble(
    ~param,          ~lower, ~upper,
    "Mu",              8000,  10000,
    "Beta",           2e-17,  5e-15,
    "BetaW",           0.85,   1.25,
    "v",                1.2,    1.5,
    "sigma",          2e-13,  4e-13,
    "phase",              6,    8.5,
    'InitInfected1' , 20000,  40000,
    'InitInfected2' ,  1500,   5000,
    'InitInfected3' ,     1,    500,
    'InitInfected4' ,     1,    500,
    'InitInfected5' ,     1,    500,
    'InitInfected6' ,     1,    500,
    'InitInfected7' ,     1,    100,
    'InitInfected8' ,  1500,   5000,
    'InitInfected9' ,     1,    500,
    'InitInfected10',     1,    500
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
  
  all_est_pars <- c(
    'Mu', 'Beta', 'BetaW', 'v', 'sigma', 'phase',
    paste0('InitInfected', 1:10), 'Delta', 'Rho', 
    'redmu', 'AlphaS'
  )
  
  all_log_params <- c(
    'Mu', 'Beta', 'BetaW', 'v', 'sigma', 
    paste0('InitInfected', 1:10), 'Delta', 'AlphaS'
  )
  
  all_logit_params <- c(
    'Rho', 'redmu'
  )
  
  est_parms <- all_est_pars[all_est_pars != prof_par]
  log_parms <- all_log_params[all_log_params != prof_par]
  logit_parms <- all_logit_params[all_logit_params != prof_par]
  
  # all_trans_pars <- est_parms[est_parms != 'phase']

  h2 <- haiti2(cutoff = 10000, measure = "log")
  start_parms <- unlist(data$starts[i, ])
  start_parms[prof_par] = prof_val

  # Get objective function to minimize
  epi_ofun <- pomp::traj_objfun(
    h2,
    est = est_parms,
    params = start_parms
  )

  log_theta <- start_parms[est_parms[est_parms %in% log_parms]]
  log_theta[log_theta == 0] <- 1
  logit_theta <- start_parms[est_parms[est_parms %in% logit_parms]]
  
  # Log-transform all parameters that need to be transformed
  theta <- c(
    log(log_theta),
    logit(logit_theta),
    start_parms[setdiff(est_parms, c(log_parms, logit_parms))]
  )[est_parms]

  # # Add "phase" parameter if needed, this variable is not transformed.
  # if (prof_par != 'phase') {
  #   theta <- c(theta, 'phase' = unlist(unname(start_parms['phase'])))
  # }

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

  all_est_pars <- c(
    'Mu', 'Beta', 'BetaW', 'v', 'sigma', 'phase',
    paste0('InitInfected', 1:10), 'Delta', 'Rho', 
    'redmu', 'AlphaS'
  )
  
  all_log_params <- c(
    'Mu', 'Beta', 'BetaW', 'v', 'sigma', 
    paste0('InitInfected', 1:10), 'Delta', 'AlphaS'
  )
  
  all_logit_params <- c(
    'Rho', 'redmu'
  )
  
  # Vector of all parameters that may be fit and transformed
  # trans_pars <- c('Mu', 'Beta', 'BetaW', 'v', 'sigma')

  # Get the names of all of the parameters that are fit for this particular
  # instance
  est_pars <- names(instance$initialization)

  # subset the transformed parameter vector using est_pars
  # trans_pars <- trans_pars[trans_pars %in% est_pars]

  # Fit the model using the given initialization
  h2_fit <- haitipkg::fit_haiti2(
    initialization = instance$initialization,
    obj_fun = instance$obj_fun,
    control = list(maxit = 20000)
  )

  estimated_vals <- h2_fit$par[est_pars]
  
  estimated_vals <- c(
    exp(h2_fit$par[all_log_params[all_log_params %in% est_pars]]),
    boot::inv.logit(h2_fit$par[all_logit_params[all_logit_params %in% est_pars]])
  )
  
  # Back transform the model parameters to natural scale
  if ('phase' %in% est_pars) {
    estimated_vals <- c(
      exp(h2_fit$par[all_log_params[all_log_params %in% est_pars]]),
      boot::inv.logit(h2_fit$par[all_logit_params[all_logit_params %in% est_pars]]),
      h2_fit$par['phase']
      )
  } else {
    estimated_vals <- c(
      exp(h2_fit$par[all_log_params[all_log_params %in% est_pars]]),
      boot::inv.logit(h2_fit$par[all_logit_params[all_logit_params %in% est_pars]])
    )
  }
  
  estimated_vals <- estimated_vals[est_pars]
  
  # Add logLik to result matrix, and rename the final vector
  results <- c('logLik' = -h2_fit$value, estimated_vals, instance$prof_val)
  names(results) <- c('logLik', est_pars, instance$prof_par)
  
  # Re-order the vector so the results are in the same order for all jobs.
  results[c(
    "logLik", "Mu", "Beta", "BetaW", "v", "sigma", "phase", 
    paste0('InitInfected', 1:10), 'Delta', 'Rho', 'redmu'
  )]
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
  account = 'stats_dept1', walltime = '60:00',
  memory = '1000m', ncpus = 1
)

resources2 <- list(
  account = 'ionides2', walltime = '60:00',
  memory = '1000m', ncpus = 1
)

submitJobs(ids = 1:350, resources = resources1)
submitJobs(ids = 351:740, resources = resources2)

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
# reg = loadRegistry('model2/profileReg')
# getStatus()
#
# Once all of the jobs are finished, we can summarize the results using:
#
# h2_profile_results = unwrap(reduceResultsDataTable())
# h2_pars = unwrap(getJobPars())
# h2_results = ijoin(h2_pars, h2_profile_results)
#
# saveRDS(
#   h2_results,
#   'model2/h2_profiles.rds')
# )
#
# We can also get information about each job, which can be helpful to pick the
# size of each run level and to estimate the cost of the RUN_LEVEL_3 search.
# For example, we can see the total / average computational time of all of the
# jobs by running:
#
# h2_stats <- unwrap(getJobStatus())
# h2_stats <- ijoin(h2_pars, h2_stats)
#
# summary(h2_stats$time.running |> as.numeric())  # unit of measurement is seconds
