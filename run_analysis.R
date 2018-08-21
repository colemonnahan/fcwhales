### This script is the central place for analyzing the FC humpback whale
### mark-recapture photoID data.

## Started 10/2017 by Cole, based off an earlier version by Noble (2015).

### Step 1: Prepare workspace and load the data for analysis
("startup.R")

## If you don't have available data you can use the simulated data to run
## the models. This will be helpful if you want to recreate the analysis
## for your own data -- start with a working version.
source("simulate_data.R")
dat <- simdata$dat
str(dat)

## Otherwise use the real data
source("load_data.R")
dat <- realdata
str(dat)

### Step 2: (____ You can skip this if already run ___). Run JAGS models,
### save model output. Fits are saved as .RDS files and loaded in
### eval_models script
source("fit_models.R")

### Step 3: Check model diagnostics
source("eval_models.R")

### Step 4: Make tables and figures for the publication
source("make_figures.R")





