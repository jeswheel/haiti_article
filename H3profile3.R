# Fitting profiles for Haiti 3: extending the profile search for parameters 
# "mu_B" and "k". 
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
#   the RUN_LEVEL == 3 version of this script took 946 compute hours. Partial
#   results can be obtained by lowering the required computational effort
#   (params: NP, NBPF, nprof, NP_EVAL, NREPS_EVAL), which can be done by
#   reducing the RUN_LEVEL to 2 or 1, or alternatively by manually changing
#   these values.
# ********
#
# This script must be ran after the completion of `H3profile.R` have been
# obtained, but can be ran simultaneously or before `H3profile2.R`. The 
# goal of this script is to extend the profile likelihood search for the 
# parameters `mu_B` and `k`. The bounds of the profile likelihood search were 
# selected by running smaller run levels. While this often results in 
# reasonable search bounds for the final profile, in this case the lower level 
# run-levels suggested that a smaller area needed to be searched for both 
# `mu_B` and `k` then what was actually needed. This script is nearly 
# identical to `H3profile2.R`, other than extending the range considered for 
# the two mentioned parameters. 

# Load Packages -----------------------------------------------------------

library(batchtools)
library(tidyverse)
library(haitipkg)

RUN_LEVEL <- 3

nprof <- switch(RUN_LEVEL, 2, 14, 20)
NBPF <- switch(RUN_LEVEL, 5, 50, 100)
NP <- switch(RUN_LEVEL, 50, 750, 1000)
NP_EVAL <- switch(RUN_LEVEL, 100, 1000, 2000)
NREPS_EVAL <- switch(RUN_LEVEL, 3, 6, 10)
SPAT_REGRESSION <- 0.05
COOLING <- 0.35

# Create Experiment Registry ----------------------------------------------

reg <- makeExperimentRegistry(
  file.dir = paste0('model3/profile_mu_B_k', RUN_LEVEL),
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

# Create Profile Design Matrix --------------------------------------------

set.seed(665544)

# Load the previous profile results, and find the maximum likelihood 
# estimate for each fixed value of the std_W profile. 
std_W_prof <- readRDS("model3/run_level_3/h3_profiles.rds") %>%
  filter(prof_var == 'std_W') %>% 
  select(logLik, std_W1) %>% 
  rename(std_W = std_W1) %>% 
  group_by(std_W) %>% 
  slice_max(order_by = logLik)

# Obtain MCAP confidence interval for std_W
mcap_res <- pomp::mcap(
  logLik = std_W_prof$logLik, 
  parameter = std_W_prof$std_W
)

# Set std_W to be the smoothed MLE from the mcap algorithm. 
std_W <- mcap_res$mle

# Create H3 model object to get fixed parameters
h3_spat <- haiti3_spatPomp()
final_pars <- coef(h3_spat)

# Create vectors for the unit and shared parameters
unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")
shared_param_names <- c(
  "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "epsilon", "k"
)

prof_params <- c("mu_B", 'k')

best_pars <- final_pars
prof_vars <- c()
for (pp in prof_params) {

  if (pp == "mu_B") {
    # Add 6 values to evaluate on the left-hand side of the profile for mu_B, 
    # using the same spacing that was used in the original profile. 
    prof_values <- seq(231.25, 325-15.625, 15.625)
  } else if (pp == 'k') {
    # Add 6 values to evaluate on the right-hand side of the profile for k, 
    # using the same spacing that was used in the original profile. 
    prof_values <- seq(150+5.416667, 150+6*5.416667, 5.416667)
  }

  if (pp %in% shared_param_names) {
    prof_cols <- matrix(rep(rep(prof_values, 10), each = nprof), ncol = 10)
    colnames(prof_cols) <- paste0(pp, 1:10)
  } else {
    prof_cols <- matrix(rep(prof_values, each = nprof), ncol = 1)
    colnames(prof_cols) <- pp
  }

  # Unit-lower bounds
  unit_lower_bounds <- c(
    'betaB1' = 2.5,
    'betaB2' = 9,
    'betaB3' = 10,
    'betaB4' = 10,
    'betaB5' = 2.1,
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
    'betaB1' = 7.5,
    'betaB2' = 30,
    'betaB3' = 40,
    'betaB4' = 40,
    'betaB5' = 9,
    'betaB6' = 45,
    'betaB7' = 15,
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
    "XthetaA" = 1e-8,
    "thetaI" = 2.5e-05,
    "lambdaR" = 0.75,
    "r" = 0.4,
    "epsilon" = 0.8,
    "k" = 30
  )

  shared_upper_bounds <- c(
    "mu_B" = 700,
    "XthetaA" = 0.15,
    "thetaI" = 1e-04,
    "lambdaR" = 2,
    "r" = 1.2,
    "epsilon" = 0.99,
    "k" = 125
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
  all_params[paste0('std_W', 1:10)] <- std_W
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  fixed_mat <- matrix(
    rep(fixed_params, nprof * length(prof_values)),
    byrow = TRUE, nrow = nprof * length(prof_values)
  )
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(h3_spat))]
  final_pars <- rbind(final_pars, guesses_all)
  prof_vars <- c(prof_vars, rep(pp, nrow(guesses_all)))
}

# Remove the first row which contains the original params as a place-holder.
final_pars <- final_pars[-1, names(coef(h3_spat))]
data_obj <- list(starts = final_pars, prof_vars = prof_vars)
# dim(final_pars)

# Define profile problem for registry -------------------------------------

# Function to sample a starting point for a given profile search
sample_starting_point <- function(data = data_obj, job, i, ...) {
  list(
    start_par = unlist(data$starts[i, ]),
    prof_var = data$prof_vars[i]
  )
}

addProblem(name = 'profile', data = data_obj, fun = sample_starting_point)

# Define algorithm for parameter estimation -------------------------------

# Function used to fit the model under the constraints of the profile search.
# If parameter "PARAM" is getting profiled over, the RW_SD is set to zero
# for "PARAM" and all other fixed parameters.
fit_model <- function(
    data, job, instance, nbpf, np, spat_regression, np_eval, nreps_eval,
    cooling, start_date = "2010-11-20", ...
) {

  # Create vectors for the unit and shared parameters
  unit_specific_names <- c("betaB", "foi_add", "aHur", "hHur")
  shared_param_names <- c(
    "mu_B", "XthetaA", "thetaI", "lambdaR", "r", "epsilon", "k"
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

pdes <- list('profile' = data.frame(i = 1:nrow(final_pars), prof_var = prof_vars))
ades <- list(
  'fitMod' = data.frame(
    nbpf = NBPF, np = NP, spat_regression = SPAT_REGRESSION,
    np_eval = NP_EVAL, nreps_eval = NREPS_EVAL, cooling = COOLING
  )
)

addExperiments(prob.designs = pdes, algo.designs = ades)

# Submit Jobs -------------------------------------------------------------

# When we ran the code, the longest running job took 4.72 hours, and the 
# average job took 3.94 hours. There are 240 total jobs. 

resources1 <- list(account = 'ionides2', walltime = '7:00:00', memory = '5000m', ncpus = 1)
submitJobs(resources = resources1)


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
# reg = loadRegistry(paste0('model3/profile_mu_B_k', RUN_LEVEL))
# getStatus()
#
# Once all of the jobs are finished, we can summarize the results using:
#
# h3_profile_results = unwrap(reduceResultsDataTable())
# h3_pars = unwrap(getJobPars())
# h3_results = ijoin(h3_pars, h3_profile_results)
#
# saveRDS(
#   h3_results,
#   paste0('model3/run_level_', RUN_LEVEL, '/', 'h3_profiles3.rds')
# )
#
# We can also get information about each job, which can be helpful to pick the
# size of each run level and to estimate the cost of the RUN_LEVEL_3 search.
# For example, we can see the total / average computational time of all of the
# jobs by running:
#
# h3_stats <- unwrap(getJobStatus())
# h3_stats <- ijoin(h3_pars, h3_stats)
#
# summary(h3_stats$time.running |> as.numeric())  # unit of measurement is seconds
