[![DOI](https://zenodo.org/badge/521646880.svg)](https://zenodo.org/badge/latestdoi/521646880)

## Haiti Article

This repository contains the necessary code to reproduce the article: "Informing policy via dynamic models: Cholera in Haiti". 

### How to recreate documents

This project contains several scripts / tasks that need to be ran in order to recreate all of the documents. These tasks should be done in the following order: 

- Run R scripts
- Compile MS (depends on results of R scripts)
- Compile SI (depends on results of MS)

Additional details breaking down each of these parts is provided below

#### Running R Scripts 

While the manuscript [`ms.Rnw`](ms.Rnw) is almost entirely self-contained, it does rely on a few data files that are created using `R` scripts.
These `R` scripts are also mostly self contained and can be ran in any order.
If a script requires running a different script first (as is the case for `H3profile` scripts), they are described as indented items below. 
Each of these scripts relies on the `R` package [`batchtools`](https://mllg.github.io/batchtools/articles/batchtools.html) in order to perform many intense calculations. 
The total computation time for these scripts is measured in tens of thousands of computing hours---which is spread over multiple machines and cores using the `batchtools` package---so be sure to be familiar with the `batchtools` package and the documentation of each of the `R` scripts before trying to run these documents. 

- [`H1profile.R`](H1profile.R): this script is used to calibrate Model~1 and to obtain confidence intervals for each of its calibrated parameters via Monte Carlo Adjusted Profile (MCAP). 
   - **Script Output**: a `batchtools` registry called `model1/profileReg_RL3`. 
   - **Target File**: `model1/run_level_3/h1_profiles.rds`.
- [`H2profile.R`](H2profile.R): this script is used to calibrate Model~2 and to obtain confidence intervals for each of its calibrated parameters via profile likelihood. 
   - **Script Output**: a `batchtools` registry called `model2/profileReg`. 
   - **Target File**: `model2/h2_profiles.rds`.
- **Model 3**: parameter estimates and confidence intervals for this model were obtained in various stages, described below. These files must be ran sequentially.
   - [`H3profile.R`](H3profile.R): this script performs MCAP for all model parameters. After running this file, it became evident that fixing the over-dispersion parameter, `std_W`, often led to better parameter estimates; this made the estimated profiles for other parameters suspect, as they often lied below the `std_W` profile. 
      - **Script Output**: a `batchtools` registry called `model3/profileReg_RL3`. 
      - **Target File**: `model3/run_level_3/h3_profiles.rds`.
   - [`H3profile2.R`](H3profile2.R): Because we found that fixing `std_W` resulted in better parameter estimates, we redo the profiles for other parameters while fixing `std_W` at it's smoothed MLE. This resulted in improved fits for each of the profiles. 
      - **Script Output**: a `batchtools` registry called `model3/profile_non_std_W_RL`. 
      - **Target File**: `model3/run_level_3/h3_profiles2.rds`.
   - [`H3profile3.R`](H3profile3.R): Finally, the bounds of the profile likelihood search were obtained by running the previous files using a smaller run-level (i.e., reducing the number of searches / number of particles used / number of IF2 iterations). These lower searches have much higher Monte Carlo error and are often not sufficient to include in final results. However, they are useful to obtain estimates of the amount of computation that is needed and determine the bounds to search over for profile confidence intervals. In this case the original bounds obtained with lower-computational searches proved to be insufficient for the parameters `mu_B` and `k`. Because of the large computational cost, this script extends the search from `H3profile2.R` for these two parameters rather than modifying the bounds in the original script and rerunning.
      - **Script Output**: a `batchtools` registry called `model3/profile_mu_B_k`. 
      - **Target File**: `model3/run_level_3/h3_profiles3.rds`.
   - [`H3_MLE_Eval.R`](H3_MLE_Eval.R): This script uses a larger number of particles and replicates to evaluate each of the parameter sets that lie within the confidence regions defined by the MCAP results above. These evaluations are done to reduce Monte Carlo variance for the final parameter estimate used in the article, while avoiding the unecessary additional computation cost to do this for parameter sets that lie in regions of low likelihood. 
      - **Script Output**: a `batchtools` registry called `model3/H3_MLE_Eval_RL`. 
      - **Target File**: `model3/run_level_3/h3_MLE_eval.rds.rds`.


Because each of these `R` scripts listed above rely on the `batchtools` package, the result of the script is a folder referred to as a `registry`. These registries need to be summarized in order to get the listed *target files*. Code to summarize the `batchtool` registries are provided at the bottom of each `R` script, but commented out since the code needs to finish executing before the summaries can be obtained. 

#### Compile MS

The manuscript must be compiled after obtaining the *target files* from the `R` scripts listed above. Once these are obtained, the manuscript is self-contained and can be created by knitting the [`ms.Rnw`](ms.Rnw) file. 

While the majority of expensive computations are done in the `R` scripts listed above, there are a few significant calculations that are done in the `ms.Rnw` file. To compile this file, it took roughly 5 hours on `36` cores, and required 14 GB of memory. 
The expensive calculations (both in terms of time-complexity and memory) in this file are wrapped in calls to the function [`pomp::bake`](https://rdrr.io/cran/pomp/man/bake.html), which caches the results so that they only need to be computed once, but are recomputed when the generating code or reproducible seed has changed. 

The file can be knitted in a few was. We ran the script on a high performance computing cluster via the command: 

```
Rscript -e "library(knitr); knit2pdf('ms.Rnw')"
```

#### Compile SI 

The manuscript ([`ms.Rnw`](ms.Rnw)) needs to be knit prior to knitting the [`si/si.Rnw`](si/si.Rnw) file. There are a few computations done in `ms.Rnw` that are later used in the supplement. 

The supplement `.Rnw` document relies on several input `.tex` files. Some of these `.tex` files are the output of a knitted `.Rnw` file. The input files are:

- [`si/inputs/header.tex`](si/inputs/header.tex)
- [`si/inputs/mod1diagram.tex`](si/inputs/mod1diagram.tex)
- [`si/inputs/mod2diagram.tex`](si/inputs/mod2diagram.tex)
- [`si/inputs/mod3diagram.tex`](si/inputs/mod3diagram.tex)
- [`si/inputs/modelDetails.tex`](si/inputs/modelDetails.tex)
- [`si/inputs/measurementModels.tex`](si/inputs/measurementModels.tex)
- [`si/initialValuesOut.tex`](si/initialValuesOut.tex)
   - This is made by knitting [`si/initialValues.Rnw`](si/initialValues.Rnw)
- [`si/calibrateMod3Out.tex`](si/calibrateMod3Out.tex)
   - This is made by knitting [`si/calibrateMod3.Rnw`](`si/calibrateMod3.Rnw`)
- [`si/ReplicateLee20Out.tex`](si/ReplicateLee20Out.tex)
   - This is made by knitting [`si/ReplicateLee20.Rnw`](`si/ReplicateLee20.Rnw`)
- [`si/confidenceIntervalsOut.tex`](si/confidenceIntervalsOut.tex)
   - This is made by knitting [`si/confidenceIntervals.Rnw`](`si/ConfidenceIntervals.Rnw`)
- [`si/inputs/paramUncertainty.tex`](si/inputs/paramUncertainty.tex)
- [`si/inputs/translateTab.tex`](si/inputs/translateTab.tex)

The [`Makefile`](si/Makefile) in the `si` folder has some commands that can knit all of the `.Rnw` files. These are relatively fast to knit, except for `si/calibrateMod3.Rnw`, which takes roughly 2.5 days on 26 cores. 

### Plos Comp Bio submission

Once the document has been knit into a `.tex` file, several steps are taken to 
prepare the document for submission to Plos Comp Bio. 
For example, no graphics should be included in the submission, and all figures
need to be Tiffs. 
Also, the `.bbl` file should be included in the `.tex` document and not be 
loaded later. 

The make file contains a command `submission` which simplifies these steps. 
Running `make submission` will create a sub-folder called `submission/` that 
contains the `ms-submission.tex` file with all of the requirements. 
It also creates the corresponding `ms-submission.pdf` so that you can check to 
make sure that everything is working as expected. 
A final step that is not currently covered by `make submission` is to convert 
the parameter table into a figure.
This is necessary because the table contains some nested tabular environments 
(created by `\multirow`). 
To do this, I simply copy `submission/ms-submission.tex`, remove fancy header /
foot and page numbers, remove the caption from the table, compile into a pdf, 
and save the single page with the parameter table to it's own `.pdf` file 
named `submission/paramTab.pdf`. 
I can then convert this to a `.tiff` by running `make submission/paramTab.tiff`. 

Finally, it is necessary to manually change the table -> figure in the 
`submission/ms-submission.tex` file, and reference the table as a figure. 

### ArXiv submission

This is primarily documented for personal benefit. The `make ArXiv` command will 
create a `.zip` file that should be able to be submitted to ArXiv. 
For this to work, you need the `.bbl` files for both `ms.tex` and `si/si.tex`, 
so these should be compiled, but not cleaned, prior to running `make ArXiv`. 
A more advanced modification to the make file would make sure these exist, but 
because this is just for personal use, it felt unecessary to spend the additional 
effort to create the make file in this case. 

### TODO: 

- [x] re-write this README with updated instructions about batchtools. 
- [x] `si/calibrateMod3.Rnw` needs to be re-factored for the new way of fitting Model 3.
- [x] Fix labels for confidence interval plots for models 2 and 3. 
- [x] Create table in confidence interval script that has all of the estimated parameters with their confidence intervals. 
- [x] Transform the confidence interval figures so that the units match the units given in the manuscript. 
- [x] Write detailed instructions on how to use to profile / batchtools scripts. 

