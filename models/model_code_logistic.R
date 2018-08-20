## Model with robust design and pr capture a logistic function of effort

## List of 6
##  $ w     : num [1:342] 1 1 1 1 1 1 1 1 1 1 ...
##  $ M     : int 342
##  $ nyrs  : int 13
##  $ nseas : int 8
##  $ x     : num [1:342, 1:13, 1:8] 0 0 0 0 0 0 0 0 0 0 ...
##  $ effort: int [1:13, 1:8] 0 0 0 0 0 0 0 0 0 1 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ year : chr [1:13] "2003" "2004" "2005" "2006" ...
##   .. ..$ month: chr [1:8] "1" "2" "3" "4" ...


model{
  effort.bar <- mean(effort)
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
        pcap[ind,yr,seas]<-a[ind,yr]*p[ind,yr,seas]
        ## Calculate probability of capture if alive
        logit(p[ind,yr,seas]) <- gamma3[seas] + gamma2*(effort[yr,seas]-effort.bar)
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

  ## Survival is constant across years and individuals
  surv ~ dunif(0,1)

  ## Priors on the logistic regression for effort
  for(seas in 1:nseas){ gamma3[seas] ~ dt(0,0.16,3) }
  gamma2~ dt(0,0.16,3)

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
