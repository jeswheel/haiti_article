# This script is used to fit Model 1.
#
# @created: Sep 2, 2022

library(haitipkg)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
library(dplyr)

cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
# cores <- 20
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)

set.seed(636813)

# Set Run Level for debug, timing, and full computation
RUN_LEVEL = 3

NP_LOCAL     <- switch(RUN_LEVEL, 100, 1000,   2000)
NMIF         <- switch(RUN_LEVEL,  20,   50,    200)
NUM_TREND    <- switch(RUN_LEVEL,   5,   10,     90)
NPROF        <- switch(RUN_LEVEL,   4,    4,     24)
NREPS_EVAL   <- switch(RUN_LEVEL,  20,   20,  cores)
NP_EVAL      <- switch(RUN_LEVEL, 200, 1000,   5000)

results <- fit_haiti1(
  NP = NP_LOCAL, NMIF = NMIF, NUM_TREND = NUM_TREND, NPROF = NPROF,
  NREPS_EVAL = NREPS_EVAL, NP_EVAL = NP_EVAL, ncores = cores
)

save_dir <- paste0("model1/run_level_", RUN_LEVEL)

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

save(results, file = paste0(save_dir, '/haiti1_fit.rda'))
