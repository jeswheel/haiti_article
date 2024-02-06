#!/bin/bash

#SBATCH --job-name=ms_RL3
#SBATCH --mail-user=jeswheel@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --account=stats_dept1
#SBATCH --partition=standard

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1

## 5GB/cpu is the basic share
#SBATCH --mem-per-cpu=1GB

## wall time hours:minutes:seconds
#SBATCH --time=05:00:00

###   Load software modules

module load R
module list

####  Commands your job should run follow this line

echo "Running on $SLURM_JOB_NODELIST"
echo "Running in $(pwd)"

Rscript -e "library(knitr); knit2pdf('ms.Rnw')"

