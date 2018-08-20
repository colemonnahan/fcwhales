## Model with robust design and pr capture a von Bertalanffy function wrt
## to effort.

### This version has a beta informative prior on survival. These parameters
## get pretty close to a confident region of (.85,.95) with mean of .929
## from Felix et al 2012 JCRM paper
## xx <- seq(.8,1, len=10000)
## alpha <- 92.9
## beta <- 10
## plot(xx, dbeta(xx, shape1=alpha, shape2=beta))
## qbeta(p=c(0.025, 0.975), shape1=alpha, shape2=beta)


model{
  ## Model dynamics
  for(ind in 1:M){ # loop over latent + observed individuals
    ## Latent state for whether whale i exists or not
    w[ind] ~ dbern(pw)
    for(yr in 1:nyrs){ # loop over years
      ## DEATH
      ad[ind,yr] ~ dbern(piad[ind,yr])
      ## BIRTH
      ab[ind,yr] ~ dbern(piab[ind,yr])
      ## ALIVE
      a[ind,yr] <- ab[ind,yr]*ad[ind,yr]*w[ind]
      for(seas in 1:nseas){ # loop over months (within year)
        ## CAPTURE
        x[ind,yr,seas] ~ dbern(pcap[ind,yr,seas])
        ## Calculate probability of capture if alive (p) or dead (0)
        pcap[ind,yr,seas] <- a[ind,yr]*p[ind,yr,seas]
        ## Calculate probability of capture if alive. This is availability
        ## times VB function relating effort to recapture probability. All
        ## three terms are between 0 and 1, so the result is a valid
        ## probability.
        p[ind,yr,seas] <- pa1[yr]*pa2[seas]*(1-exp(-k*effort[yr,seas]))
        ## Posterior predictive check: bernoulli draw for whether whales alive
        ## were seen with probability of p across seasons and
        ## years. Multiplying by a zeroes it out if animal is dead in that
        ## year
        post_predict_all[ind,yr,seas] ~ dbern(a[ind,yr]*p[ind,yr,seas])
      } # end loop months
    } # end loop years

    ## Initiate prob for birth (piab) and prob of death (piad)
    piab[ind,1] <- zeta0
    piad[ind,1] <- 1
    ## sazero tracks whether available to be born (1) or not (0), the
    ## \psi_{i,j} in the paper. I.e., whether was recruited previously.
    sazero[ind,1] <- 1
    for(yr in 2:nyrs){ # loop over year
      sazero[ind,yr] <- sazero[ind,yr-1]*(1-ab[ind,yr-1])
      ## probability of being born in this period (=1 if already born)
      piab[ind,yr] <-
        ## zeta needs to be 1 in the very last year to ensure that it was
        ## born before this period, hence the if statement.
        ifelse(yr==nyrs,
               sazero[ind,yr]*1.0 + (1-sazero[ind,yr]),
               sazero[ind,yr]*zeta + (1-sazero[ind,yr]))
      piad[ind,yr] <-
        ad[ind,yr-1]*(ab[ind,yr-1]*surv + (1-ab[ind,yr-1]))
      ## piad = probability of surviving time period yr
      ## = (alive previous year) * (recruited + survived last year OR not recruited last year)
    }
  }  # end loop over individuals
  ## End of model dynamics

  ## ----- Parameter priors  -----
  ## We estimate two probabilities of birth. The first year is different
  ## b/c it reprents all births before the study (it will thus be
  ## higher). The second is a fixed effects for the rest, which represents
  ## the probability of birth in a year. We put uniform on both.
  zeta0 ~ dunif(0,1 ) # initial
  zeta ~ dunif(0,1) # annual

  ## Survival is constant across years and individuals.
  surv ~ dbeta(92.9, 10)

  ## Priors on the VB relationship with effort
  ## Seasonal fixed effects
  for(seas in 1:nseas){ pa2[seas] ~ dbeta(1,1) }
  ## Random annual effects. Hyperpriors which give about a uniform on the
  ## inverse logit (probability) scale. I.e., resemble U(0,1).
  mu_avail ~ dnorm(0,1) #dnorm(2.2,1/(.05*.05))
  sd_avail ~ dnorm(0,2)T(0,) #dunif(0,5)
  prec_avail <- 1/(sd_avail*sd_avail)
  ## Likelihood of annual availability random effects
  for(yr in 1:nyrs){
    tau_avail[yr] ~ dnorm(mu_avail, prec_avail)
    logit(pa1[yr]) <- tau_avail[yr]
  }
  ## Fixed effect k is the growth rate. As goes to zero it is linear, if it
  ## is big it saturates quickly.
  k~dunif(0,.5) #dnorm(0.1,1/(0.1*0.1))T(0,)

  ## Uniform prior on probability of latent whales
  pw ~ dbeta(1,1)

  ## Calculate derived quantities
  for(yr in 1:nyrs){
    ## Annual abundance
    N[yr]<-sum(a[1:M,yr])
  }
  ## Posterior predictive checks. Note that post_predict_all is generated
  ## above for convenience.
  for(yr in 1:nyrs){
    ## Collapse the posterior predictive distribution into # unique whales
    ## seen per year
    for(ind in 1:M){
      ## collapse to whether ind observed in at least one season in that year
      post_predict2[ind,yr] <- ifelse(sum(post_predict_all[ind,yr,])>0,1,0)
    }
    ## sum over total unique individuals seen in at least one season
    post_predict[yr] <- sum(post_predict2[,yr])
  }
  TotalN <- sum(w[1:M])
}# end model
