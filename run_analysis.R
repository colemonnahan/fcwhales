### This script is the central place for analyzing the FC humpback whale
### mark-recapture photoID data.

## Started 10/2017 by Cole, based off an earlier version by Noble (2015).

### Step 1: Load data and prepare workspace
source("startup.R")

### Step 2: Run JAGS models, save model output. Only run these if needed,
### as they are saved as .RDS files and read back in in load_results.R
## This old file does some simulation testing
source("run_models_simulated.R")
## Otherwise this one
source("run_models.R")

### Step 3: Create exploratory plots and check model diagnostics
source("load_results.R")
source("make_plots.R")
source("eval_models.R")

### Step 4: Make tables and figures for the publication
source("make_figures.R")





