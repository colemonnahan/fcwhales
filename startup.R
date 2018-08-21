library(ggplot2)
library(plyr)
library(reshape2)
library(rjags)
library(coda)
library(boot)
library(shinystan)
library(rstan)
library(R2jags)

message("Clearing workspace..")
rm(list=ls())
message("Loading libraries, preparing workspace, loading data..")

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
vb.inits <- function(dat, tvzeta=FALSE){
  f <- function(){
    x <- list(surv=runif(1,0,1),
              mu_avail=rnorm(1, 0,2),
              sd_avail=runif(1,0,2),
              tau_avail=rnorm(dat$nyrs,0,1),
              k=runif(1, 0, .5),
              zeta0=runif(1,0,1),
              ad=matrix(nrow=dat$M, ncol=dat$nyrs, data=1),
              ab=matrix(nrow=dat$M, ncol=dat$nyrs, data=1))
    ## Two cases for zeta which change the parameters in JAGs. It can be
    ## constant:
    if(!tvzeta) {
      x$zeta <- runif(1,0,1)
    } else {
      ## or time varying (ie random effects)
      x$tau_zeta <- runif(dat$nyrs, -1,1)
      x$mu_zeta <- runif(1, -1,1)
      x$sd_zeta <- runif(1, 0,1)
    }
    return(x)
  }
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
run_vb <- function(datlist, n.iter, thin, n.chains, surv='uniform',
                   tvzeta=FALSE){
  ## these are necessary to get jags to run in parallel
  list2env(datlist, envir=globalenv() )
  list2env(list(n.iter=n.iter), envir=globalenv() )
  list2env(list(thin=thin), envir=globalenv() )
  surv <- match.arg(surv, c('uniform', 'informative'))
  if(surv=='uniform'){
    ## uniform survival and constant zeta
    model.file <- 'models/model_code_vb.R'
  } else if(!tvzeta){
    ## informative survival but constant zeta
    model.file <- 'models/model_code_vb_surv.R'
  } else {
    ## informative survival and random effects on zeta
    model.file <- 'models/model_code_vb_surv_tvzeta.R'
  }
  pars <-
    c('surv', 'mu_avail', 'sd_avail', 'pa1', 'k',
      'pa2', 'pw', 'zeta0', 'zeta', 'N',
      'deviance', 'TotalN', 'post_predict' )
  if(tvzeta) pars <- c(pars, 'tau_zeta', 'sd_zeta', 'mu_zeta')
  mod.vb.fit <-
    jags.parallel(model.file=model.file, data=datlist,
                  inits=vb.inits(datlist, tvzeta=tvzeta),
                  n.iter=n.iter,
                  n.thin=thin,
                  n.chains=n.chains,
                  ##n.burnin=n.adapt,
                  parameters.to.save=pars)
  return(mod.vb.fit)
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

### These functions are used to simulate data for testing the model
#' Calculate a probability of recapture based on effort and a parameters.
#'
#' Uses a von Bertalanffy growth formula a*(1-exp(-b*effort)) to predict
#'   probability of capture given effort. Constraining 0<a<1 and b>0
#'   ensures the return value are in [0,1]. Also note that this function
#'   can be flat (b goes to Inf) or linear (b goes to 0).
#'
#' @param e Effort vector (days of effort in a month)
#' @param a Equivalent to Linf parameter of the VBLG formula (saturation
#'   point), in interval [0,1].
#' @param b Equivalent to k parameter (controls how fast saturates)
#' @return A vector of probablity of captures
effort.fn <- function(e, pa1, pa2, k){
  stopifnot(0 < pa1 & pa1 <= 1)
  stopifnot(0 < pa2 & pa2 <= 1)
  pcap <- pa1*pa2*(1-exp(-k*e))
  return(pcap)
}
## ## Quick code to check this function is working
## e <- seq(0, 20, len=1000)
## k <- .1
## xx <- ldply(c(1, .7, .9, 1), function(pa1)
##   ldply(c(.1, .2, .5, .15), function(pa2)
##     data.frame(pa1=pa1, pa2=pa2, e=e, pcap=effort.fn(e,pa1,pa2,k), group=paste(pa1,pa2, sep="_"))))
## ggplot(xx, aes(e, pcap, group=group, color=factor(pa2))) + geom_line() +
##   facet_wrap("pa1") + ylim(0,1)

#' Simulate a data set that generally mimics the FC humpback whale
#' and data collection scheme.
#'
#' @param N The true population size (number of whales)
#' @param par A list of true parameters
#' @param effort A matrix (years x months) of effort (number of days) in a
#'   month.
#' @param nyrs The number of years
#' @param nseas The number of seasons (months) within a year
#' @param seed Random number seed to set, if given
#'
#' @return A list containing data processed for use in the JAGS
#'   model. I.e., it has the same structure as the real data.
simulate_data <- function(N, pars, effort, nyrs=13, nseas=8, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  stopifnot(nyrs == nrow(effort))
  stopifnot(nseas == ncol(effort))
  ## Kinda tricky to simulate the process part of the model, where I want
  ## to specify the number of whales alive at some point in the
  ## study. There are no latent states in reality!
  piab <- pars$piab
  surv <- pars$surv
  pa1 <- pars$pa1
  pa2 <- pars$pa2
  k <- pars$k
  stopifnot(length(pa1)==nyrs)
  stopifnot(length(pa2)==nseas)
  ## The empty capture/recpature matrix and probabliity of captures array
  x.true <- array(NA, dim= c(N, nyrs, nseas))
  pcap <- matrix(NA, nyrs, nseas)
  alive <- matrix(0, N, nyrs) # matrix of whether alive
  probs.birth <- c(piab*20, rep(piab, nyrs-1))
  probs.birth <- probs.birth/sum(probs.birth)
  ## Loop through each individual whale and generate when it was alive
  ## (stochastic birth and death) and from that stochastically sample
  ## capture and recaptures.
  for(ind in 1:N){
    ## First probability here is being born in the last 20 years or on the
    ## year of the first sampling period. This ensures most whales were born
    ## before sampling, but some can still be born during the nyrs time
    ## period.
    year.born <- which(as.numeric(rmultinom(n=1, size=1, prob=probs.birth))==1)
    ## After being recruited they can die, so determine if that happens
    ## withing sampling period
    test <- rbinom(n=(nyrs-year.born), size=1, prob=surv)
    year.died <- year.born+which(test==0)[1] # NA if survived whole period
    ## Years available to be observed (= ab*ad)
    alive[ind, year.born:min(year.died, nyrs, na.rm=TRUE)] <- 1
  }
  ## I split this apart from process otherwise the sampling affected the
  ## population dynamics due to random seeds. This way is more easily
  ## reproducible
  for(ind in 1:N){
    ## Now the observation part of the model.
    for(yr in 1:nyrs){
      for(seas in 1:nseas){
        pcap[yr,seas] <-
          effort.fn(e=effort[yr,seas], pa1=pa1[yr], pa2=pa2[seas], k=k)
        ## Try to capture the whale in yr/seas
        x.true[ind,yr,seas] <- rbinom(n=1, size=1, prob=alive[ind,yr]*pcap[yr,seas])
      }
    }
  }
  abundance <- apply(alive, 2, sum) # abundance by year
  observed <- apply(apply(x.true, c(1,2) , max), 2, sum) # unique whales by year
  year.first.obs <- apply(apply(x.true, c(1,2) , max), 1, function(x)
    which(x==1)[1])
  observed.cumsum <- as.numeric(cumsum(table(factor(year.first.obs, levels=1:nyrs))))
  truth <- list(abundance=abundance, pcap=pcap, observed=observed,
  observed.cumsum=observed.cumsum, pars=pars)
  ## Prepare data for output
  which.unobserved <- which(apply(x.true, 1, sum)==0)
  if(length(which.unobserved)>0){
    x.true <- x.true[-which.unobserved,,]
    print(paste(length(which.unobserved), "whales unobserved"))
  }
  M <- floor(1.5*N) # number of allowable latent whales
  ## Build new x matrix mimicing how the analyst would do it
  x <- array(0, dim= c(M, nyrs, nseas))
  x[1:nrow(x.true),,] <- x.true
  w <- c(rep(1, N), rep(NA, (M-N))) # latent whales to initialize
  dat <- list('w'=w, 'M'=M, 'nyrs'=nyrs, 'nseas'=nseas, 'x'= x,
              'effort'=effort)
  return(list(dat=dat, truth=truth))
}

## quick function to explore simulated data properties
plot.simdat <- function(sim){
  n <- sim$dat$nyrs
  plot(0, type='n', xlim=range(1:n), ylim=c(0,150), xlab="Year", ylab='Whales')
  lines(1:n, sim$truth$abundance, lwd=3)
  points(1:n, sim$truth$observed, lty=1, pch=16, col=4, type='b')
  legend('topleft', legend=c("Sim. Abundance", "Sim. whales obs."),
         lty=c(1,2), pch=c(NA, 16), col=c(1,4),
         lwd=c(2,1))
}
