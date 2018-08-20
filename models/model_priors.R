model{
  x~dt(0,1/(20^2),1)
  for(seas in 1:8){ pa2[seas] ~ dbeta(1,1) }
  ## Random annual effects
  mu_avail ~ dnorm(2,1/(.1*.1))
  sd_avail ~ dnorm(0,10)T(0,) # half-normal
  prec_avail <- 1/(sd_avail*sd_avail)
  ## Likelihood of annual availability random effects
  for(yr in 1:20){
    tau_avail[yr] ~ dnorm(mu_avail, prec_avail)
    logit(pa1[yr]) <- tau_avail[yr]
  }
  ## Fixed effect k is the growth rate. As goes to zero it is linear, if it
  ## is big it saturates quickly.
  k~dnorm(0,1)T(0,)
  }
