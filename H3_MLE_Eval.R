# Finding and evaluating the MLE for model 3 using the profile likelihood 
# searches. 
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
#   the RUN_LEVEL == 3 version of this script took 320 compute hours. Partial
#   results can be obtained by lowering the required computational effort
#   (params: NP, NBPF, nprof, NP_EVAL, NREPS_EVAL), which can be done by
#   reducing the RUN_LEVEL to 2 or 1, or alternatively by manually changing
#   these values.
# ********
#
# In this script we use the results of the profile likelihood searches to find 
# and evaluate the MLE of Model 3. The profile likelihood searches (done using 
# the scripts H3profile.R, H3profile2.R, and H3profile3.R) contain the results 
# of a large number of searches and computation. These scripts result in 
# data.frame objects, where each row is a different parameter estimate. 
# Here, we consider only parameter estimates that are within the smoothed MCAP
# confidence intervals for all profile likelihoods, and evaluate their 
# likelihood. The time-complexity of each job is roughly linear in the 
# parameters NP_EVAL, and NREPS_EVAL. 


# Load Packages -----------------------------------------------------------

library(batchtools)
library(tidyverse)
library(haitipkg)

RUN_LEVEL <- 3

NP_EVAL <- switch(RUN_LEVEL, 100, 1000, 5000)
NREPS_EVAL <- switch(RUN_LEVEL, 3, 6, 15)

# Create Experiment Registry ----------------------------------------------

reg <- makeExperimentRegistry(
  file.dir = paste0('model3/H3_MLE_Eval_RL', RUN_LEVEL),
  seed = 739164,
  packages = c("spatPomp", 'haitipkg', 'pomp')
)

# Create Profile Design Matrix --------------------------------------------

# Load the first profiles, primarily used for getting MLE of std_W
H3_profiles1 <- readRDS('model3/run_level_3/h3_profiles.rds') %>% 
  filter(prof_var == 'std_W')

# Load the second profiles, which have the majority of the results 
H3_profiles2 <- readRDS("model3/run_level_3/h3_profiles2.rds")

# Profiles3 contain additional searches that were necessary for getting a 
# complete profile for parameters "mu_B" and "k", which was done by extending 
# the range of values considered for these parameters in their respective 
# profiles. 
H3_profiles3 <- readRDS("model3/run_level_3/h3_profiles3.rds")

# Combine all results into a single data.frame object. 
H3_profiles <- bind_rows(H3_profiles1, H3_profiles2, H3_profiles3)

# Get the names of all of the profile parameters 
prof_params <- unique(H3_profiles$prof_var)


# Get 95% confidence intervals --------------------------------------------
# In this section, we compute the 95% confidence intervals for each of the 
# profiled parameters. This is used so that we evaluate all parameter 
# estimates that lie in each interval. 

# Convert profile data.frame to "long" format, and select only the parameter 
# set that corresponds to the maximum likelihood for each value of each 
# profiled parameter. 
H3_profiles_long <- H3_profiles %>%
  select(prof_var, logLik, logLik_se, all_of(paste0(prof_params, '1'))) %>%
  rename_with(~ sub("1$", "", .x)) %>%
  select(-epsilon, -XthetaA) %>% 
  pivot_longer(
    cols = -c(prof_var, logLik, logLik_se),
    names_to = 'variable',
    values_to = 'value'
  ) %>%
  filter(variable == prof_var) %>%
  group_by(prof_var, value) %>%
  slice_max(order_by = logLik, n = 1)

# Do the same for XthetaA and epsilon, but they require a slightly different 
# treatment because the profile was computed on a transformed scale, due to 
# the MLE lying near the boundary of possible values (i.e., XthetaA >= 0, and 
# 0< epsilon <= 1). 
H3_remaining_long <- H3_profiles %>% 
  select(prof_var, logLik, logLik_se, epsilon1, XthetaA1) %>% 
  rename(epsilon = epsilon1, XthetaA = XthetaA1) %>% 
  mutate(XthetaA = sqrt(XthetaA)) %>%
  mutate(epsilon = epsilon^2) %>%
  pivot_longer(
    cols = -c(prof_var, logLik, logLik_se),
    names_to = 'variable',
    values_to = 'value'
  ) %>%
  filter(variable == prof_var) %>%
  group_by(prof_var, value) %>%
  slice_max(order_by = logLik, n = 2)

# Combined both long-versions of the results into a single data.frame
H3_profiles_long <- bind_rows(
  H3_profiles_long, H3_remaining_long
)

# Loop through all of the profile parameters and compute the confidence 
# interval. 
all_ci <- list()
for (p in prof_params) {
  
  mcap_tmp <- mcap(
    logLik = H3_profiles_long %>% filter(prof_var == p) %>% pull(logLik),
    parameter = H3_profiles_long %>% filter(prof_var == p) %>% pull(value),
    span = case_when(
      p == 'r' ~ 1,
      p == 'XthetaA' ~ 0.7,
      TRUE ~ 0.75
    )
  )
  
  tmp_ci <- mcap_tmp$ci
  names(tmp_ci) <- c("lower", "upper")
  
  all_ci[[p]] <- c(tmp_ci, 'mle' = mcap_tmp$mle)
}

# Save the confidence intervals as a data.frame object
all_ci <- all_ci |> as.data.frame() |> t() |> as.data.frame()
all_ci$prof_var <- rownames(all_ci)

# Transform epsilon and XthetaA to natural scale
all_ci <- all_ci %>% 
  mutate(
    lower = case_when(
      prof_var == 'epsilon' ~ sqrt(lower), 
      prof_var == 'XthetaA' ~ lower ^ 2,
      TRUE ~ lower
    ),
    upper = case_when(
      prof_var == 'epsilon' ~ sqrt(upper), 
      prof_var == 'XthetaA' ~ upper ^ 2,
      TRUE ~ upper
    ),
    mle = case_when(
      prof_var == 'epsilon' ~ sqrt(mle), 
      prof_var == 'XthetaA' ~ upper ^ 2,
      TRUE ~ mle
    )
  )

# Find all parameters that lie within the 95% confidence interval. 
eval_set <- H3_profiles %>% 
  filter(
    std_W1 < all_ci['std_W', 'upper'],
    std_W1 > all_ci['std_W', 'lower'],
    mu_B1 < all_ci['mu_B', 'upper'],
    mu_B1 > all_ci['mu_B', 'lower'],
    r1 < all_ci['r', 'upper'],
    r1 > all_ci['r', 'lower'],
    thetaI1 > all_ci['thetaI', 'lower'],
    thetaI1 < all_ci['thetaI', 'upper'],
    epsilon1 < all_ci['epsilon', 'upper'],
    epsilon1 > all_ci['epsilon', 'lower'],
    XthetaA1 > all_ci['XthetaA', 'lower'],
    XthetaA1 < all_ci['XthetaA', 'upper'],
    lambdaR1 < all_ci['lambdaR', 'upper'],
    lambdaR1 > all_ci['lambdaR', 'lower'],
    k1 > all_ci['k', 'lower'],
    k1 < all_ci['k', 'upper']
  ) %>% 
  select(-job.id, -problem, -algorithm, -i, -nbpf, -np, -spat_regression, 
         -cooling, -logLik, - logLik_se, -np_eval, -nreps_eval)

# Define profile problem for registry -------------------------------------

# Function to sample a starting point for a given profile search
sample_starting_point <- function(data = eval_set, job, i, ...) {
  list(
    start_par = data %>% 
      slice(i) %>% 
      select(-prof_var) %>% 
      unlist(),
    prof_var = data %>% 
      slice(i) %>% 
      pull(prof_var)
  )
}

addProblem(name = 'eval', data = eval_set, fun = sample_starting_point)

# Define algorithm for parameter estimation -------------------------------

# Function used to evaluate the model given the input dataset. The model is 
# evaluated using the block particle filter.
eval_model <- function(
    data, job, instance, np_eval, nreps_eval, start_date = "2010-11-20", ...
) {
  
  # Create the model 
  h3_spat <- haiti3_spatPomp(start_date = start_date)
  
  # Set model coefficients to input parameter value
  coef(h3_spat) <- instance$start_par
  
  # Evaluate using NREPS_EVAL replicates of the block particle filter, with 
  # NP_EVAL particles for each replicate. 
  evals <- replicate(
    nreps_eval, logLik(bpfilter(h3_spat, Np = np_eval, block_size = 1))
  )
  
  # estimate the log-likelihood and corresponding SE 
  ll <- pomp::logmeanexp(evals, se = TRUE)
  names(ll) <- c("logLik", "logLik_se")
  
  # Return log-likelihood estimate with parameter value. 
  c(ll, instance$start_par)
}

addAlgorithm(name = 'evalMod', fun = eval_model)

# Completing the experiment -----------------------------------------------

pdes <- list('eval' = data.frame(i = 1:nrow(eval_set), prof_var = eval_set$prof_var))
ades <- list(
  'evalMod' = data.frame(
    np_eval = NP_EVAL, nreps_eval = NREPS_EVAL
  )
)

addExperiments(prob.designs = pdes, algo.designs = ades)

# Submit Jobs -------------------------------------------------------------

# The maximum run-time of the 168 jobs was 2.20 hours, with a mean time of 
# 1.91 hours

resources1 <- list(account = 'stats_dept1', walltime = '3:00:00', memory = '5000m', ncpus = 1)
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
# reg = loadRegistry(paste0('model3/H3_MLE_Eval_RL', RUN_LEVEL))
# getStatus()
#
# Once all of the jobs are finished, we can summarize the results using:
#
# h3_eval_results = unwrap(reduceResultsDataTable())
# h3_pars = unwrap(getJobPars())
# h3_results = ijoin(h3_pars, h3_eval_results)
#
# saveRDS(
#   h3_results,
#   paste0('model3/run_level_', RUN_LEVEL, '/', 'h3_MLE_eval.rds')
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
