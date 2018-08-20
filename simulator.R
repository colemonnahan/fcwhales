### This file contains functions to simulate data for testing the model

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
#'
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
  observed <- apply(apply(dat$x, c(1,2) , max), 2, sum) # unique whales by year
  plot(0, type='n', xlim=range(1:n), ylim=c(0,150), xlab="Year", ylab='Whales')
  lines(1:n, sim$truth$abundance, lwd=3)
  points(1:n, sim$truth$observed, lty=1, pch=16, col=4, type='b')
  ## lines(1:n, sim$truth$observed.cumsum, lty=2)
  lines(1:length(observed), observed, pch=15, type='b', col='red')
  legend('topleft', legend=c("Sim. Abundance", "Sim. whales obs.",
                             "Real whales obs."),
         lty=c(1,2,2), pch=c(NA, 16, 16), col=c(1,4,2),
         lwd=c(2, 1,1))
}
