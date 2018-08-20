library(ggplot2)
library(plyr)
library(reshape2)
library(rjags)
load.module("dic")  # for monitoring deviance
library(coda)
library(boot)
library(shinystan)
library(rstan)
library(R2jags)

message("Clearing workspace..")
rm(list=ls())
message("Loading libraries, preparing workspace, loading data..")
source("simulator.R")

## function to convert coda output so it works with rstan (list to array)
convert_array <- function(x){
  y <- array(NA, dim=c(dim(x[[1]])[1], length(x), dim(x[[1]])[2]))
  dimnames(y) <- list(NULL, NULL, unlist(dimnames(x[[1]])[2]))
  for(i in 1:length(x)) y[,i,] <- x[[i]]
  y
}
## Plotting settings
ggwidth <- 7 # inches
ggheight <- 5

months <- c('Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun')# 1:8 # labels used throughout
## This is the calendar year of December of the effort year. So the
## 2011/2012 season is coded as 2011.
years <- 2002:2016
nyrs <- length(years)
nseas <- length(months)
n.timesteps <- nyrs*nseas

### Read in raw data from files
message("Loading sightings data...")
captures <- read.csv("data/capture_histories_2017_11_28.csv") # capture histories
n.obs <- dim(captures)[1] # number of unique whales
## Make sure each year has 8 columns and each month has 15 years
tt <- names(captures)
year.check <- any(nseas !=as.numeric(table(as.numeric(substr(tt,2,5)))))
if(year.check) stop("Wrong number of seasons for years in capture data")
month.check <- any(nyrs != as.numeric(table(substr(tt,7,9)) ))
if(month.check) stop("Wrong number of years in a season in capture data")
if(n.timesteps != dim(captures)[2])
  stop("wrong number of columns in capture history data file")

### Process data for modeling and plotting
## Add extra rows to make space for the latent whales
M <- 250#2*n.obs+50
post <- matrix(nrow=M, ncol=n.timesteps, data=0)
post[1:n.obs, 1:n.timesteps] <- as.matrix(captures)
## Convert capture histories into array form
x <- array(data=NA, dim=c(M, nyrs, nseas))
dimnames(x) <- list("ID"=1:dim(x)[1], "year"=years, "month"=months)
for(i in 0:(nyrs-1) ){
	x[ ,i+1 , ] <- post[, (i*nseas + 1):(i*8 + nseas)]
}
if(sum(is.na(x)) > 0)
  stop("Issue converting captures to array -- NAs present")

message("Loading effort data...")
## Effort:year x month. Each cells is the number of days of effort.
pre.effort <- read.table("data/effort_2017_11_28.csv", header=F, sep=',')
## transpose so month x year
effort <- as.matrix(pre.effort)
dimnames(effort) <- list("year"=years, "month"=months)
effort.long <- melt(effort, value.name='effort')


message("Checking for data issues...")
captures.long <- melt(x)
captures.long <- ddply(captures.long, .(month, year), summarize, unique.whales=sum(value))
## Use this to check back with original data if needed
## captures.long <- captures.long[order(captures.long$year, captures.long$month),]
obs.effort <- merge(effort.long, captures.long, by=c("year", "month"))
effort.check <- droplevels(subset(obs.effort, effort == 0 & unique.whales > 0))
if(nrow(effort.check)>0){
  print(effort.check)
  warning("These months have 0 effort but observed whales")
}

## This is whether each whale was seen in a year
seen.years <- apply(x, c(1,2), max)
## Now get the first year seen
temp <- as.numeric(apply(seen.years, 1, function(ii) which(ii==1)[1]))
## Drop the NAs and convert to year.
first.year.seen <- years[temp[1:nrow(captures)]]
## Trick to get table to count 0 whales when year is missing
first.year.seen <- factor(first.year.seen, levels=years)
x1 <- as.numeric(table(first.year.seen))
x2 <- as.numeric(apply(apply(x, c(1,2) , max), 2, sum)) # unique whales by year
e <- as.numeric(apply(effort, 1, sum))
whales.table <- data.frame(year=years+1, effort=e, unique.whales=x2, new.whales=x1)

### No longer using this.
## ## start months and end months for each visit
## non.visit <- read.table("data/not_visited15.csv", header=T, row.names=1, sep=',')
## if(nseas != dim(non.visit)[1])
##   stop("wrong number of rows in non_visit data file")
## NV <- as.matrix(non.visit)
## dimnames(NV) <- list("month"=months, "year"=years)

## Prep the latent whale JAGS inputs
w <- rep(NA, length=M)
w[1:n.obs] <- 1
## Real data for all models
dat <- list('w'=w, 'M'=M, 'nyrs'=nyrs, 'nseas'=nseas, 'x'= x,
                 'effort'=effort)

message("Loading initialization functions...")
## Converts MCMC output list into dataframe with better column names for
## plotting.
convert_to_df <- function(fit, name=NULL){
  fit <- as.mcmc(fit)
  post <- data.frame(do.call(rbind, fit))
  new <- gsub('[.]$', '', names(post))
  new <- gsub('[.]','_', new)
  names(post) <- new
  if(!is.null(name)) post$model <- name
  return(post)
}


## Init functions.
vb.inits <- function(dat){
  f <- function()
  list(surv=runif(1,0,1),
       mu_avail=rnorm(1, 0,2),
       sd_avail=runif(1,0,2),
       tau_avail=rnorm(dat$nyrs,0,1),
       k=runif(1, 0, .5),
       zeta0=runif(1,0,1),
       zeta=runif(1,0,1),
       ad=matrix(nrow=dat$M, ncol=dat$nyrs, data=1),
       ab=matrix(nrow=dat$M, ncol=dat$nyrs, data=1))
  return(f)
}
logistic.inits <- function(dat){
  f <- function()
    list(surv=runif(1,0,1),
         zeta0=runif(1,0,1),
         zeta=runif(1,0,1),
         gamma3=runif(dat$nseas,-2,2),
         gamma2=runif(1, -2, 2),
         ad=matrix(nrow=dat$M, ncol=dat$nyrs, data=1),
         ab=matrix(nrow=dat$M, ncol=dat$nyrs, data=1))
}

## quick function to run the logistic model
run_logistic <- function(datlist, n.iter, thin, n.chains){
  ## these are necessary to get jags to run in parallel
  list2env(datlist, envir=globalenv() )
  list2env(list(n.iter=n.iter), envir=globalenv() )
  list2env(list(thin=thin), envir=globalenv() )
  mod.logistic.fit <-
    jags.parallel(model.file='models/model_code_logistic.R',
                  data=datlist,
                  n.iter=n.iter,
                  n.thin=thin,
               inits=logistic.inits(datlist),
               n.chains=n.chains,
               #n.burnin=n.adapt,
               parameters.to.save=c('surv', 'gamma2', 'gamma3',
                  'zeta0', 'zeta', 'N', 'pw', 'deviance' ,'TotalN', 'post_predict'))
  return((mod.logistic.fit))
}


## quick function to run the vb model
run_vb <- function(datlist, n.iter, thin, n.chains, surv='uniform'){
  ## these are necessary to get jags to run in parallel
  list2env(datlist, envir=globalenv() )
  list2env(list(n.iter=n.iter), envir=globalenv() )
  list2env(list(thin=thin), envir=globalenv() )
  surv <- match.arg(surv, c('uniform', 'informative'))
  model.file <- ifelse(surv=='uniform', 'models/model_code_vb.R',
                       'models/model_code_vb_surv.R')
  mod.vb.fit <-
    jags.parallel(model.file=model.file, data=datlist,
                  inits=vb.inits(datlist),
                  n.iter=n.iter,
                  n.thin=thin,
               n.chains=n.chains,
               #n.burnin=n.adapt,
               parameters.to.save=c('surv', 'mu_avail', 'sd_avail', 'pa1', 'k',
                                    'pa2', 'pw', 'zeta0', 'zeta', 'N',
                                    'deviance', 'TotalN', 'post_predict' ))
  return((mod.vb.fit))
}


add.polygon <- function(x, y, z, alpha.level, alpha.min=0, alpha.max=1){
    ## Pass this a series of years (x) and a matrix of trajectories (df) and it
    ## adds a polygon to the current plot with whose area contains 1-alpha.level
    library(reshape2)
    alpha.level <- sort(alpha.level)
    for(i in 1:length(alpha.level)){
        alpha <- alpha.level[i]
        alpha.col <- alpha.min+alpha*(alpha.max-alpha.min)
        col.poly <- rgb(1-alpha.col,1-alpha.col,1-alpha.col, alpha=1)
        quantiles.temp <-  as.matrix(t(apply(z, 2, quantile,
                                             probs=c(alpha/2,1-alpha/2),name=F, na.rm=T)))
        polygon(x=c(x, rev(x)), y=c(quantiles.temp[,1], rev(quantiles.temp[,2])),
                col=col.poly, border=NA)
    }
}
