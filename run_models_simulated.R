## This file generates two simulated data sets and fits models to them,
## writing all to the results folder.

## The goal is to verify the models are performing as expected under known
## circumstances.
message("Simulating data sets...")
set.seed(2131424)

## Simulated dataset 1: Bias check. Throw a ton of data at it and make sure
## model recovers parameters. Effort is highly variable to try and capture
## shape of effort curve.
nyrs2 <- 30; nseas2 <- 8
effort.sim1 <- matrix(data=sample(x=0:20, size=nyrs2*nseas2, replace=TRUE) , nrow=nyrs2, nseas2)
## draw  centered on the prior
pa1 <- rnorm(nyrs2, 2,.3)
pa1[5] <- .5 # one really low year to see if it can recover this
pars <- list(piab=0.2, surv=.9, pa1=inv.logit(pa1),
             pa2=c(.05,.1, .2, .4, .5, .3, .2, .08),
             k=.2)
sim1 <- simulate_data(N=250, pars=pars, effort=effort.sim1, seed=152,
                      nyrs=nyrs2, nseas=nseas2)
sim1$dat$M <- 275
saveRDS(sim1, 'results/sim1.RDS')
png("plots/simdata1.png", width=7, heigh=5, units='in', res=500)
plot.simdat(sim1)
mtext("Data rich simulated data properties", line=1.5, cex=1.5)
dev.off()

## Simulated dataset 2: Mimic real data. Pass real effort matrix. Make sure
## can recover parameters
nyrs2 <- nyrs; nseas2 <- nseas
## Use the real effort data
pa1 <- rnorm(nyrs2, 2,.3)
pa1[5] <- .5 # one really low year to see if it can recover this
effort.sim2 <- dat$effort
pars <- list(piab=0.2, surv=.9, pa1=inv.logit(pa1),
             pa2=c(.05,.1, .2, .4, .5, .3, .2, .08),
             k=.2)
sim2 <- simulate_data(N=150, pars=pars, effort=effort.sim2, seed=152,
                      nyrs=nyrs2, nseas=nseas2)
saveRDS(sim2, 'results/sim2.RDS')
png("plots/simdata2.png", width=7, heigh=5, units='in', res=500)
plot.simdat(sim2)
mtext("Realistic simulated data properties", line=1.5, cex=1.5)
dev.off()


## Fit Noble's original logistic regression model to both
th <- 1; ni <- 10*th; nc <- 2; (ni/2/th*nc)
th <- 50; ni <- 200*th; nc <- 10; (ni/2/th*nc)

start <- Sys.time()
fit.sim1.logistic <-
  run_logistic(dat=sim1$dat, n.iter=ni, thin=th, n.chains=nc)
saveRDS(fit.sim1.logistic, file='results/fit.sim1.logistic.RDS')

fit.sim2.logistic <-
  run_logistic(dat=sim2$dat, n.iter=ni, thin=th, n.chains=nc)
saveRDS(fit.sim2.logistic, file='results/fit.sim2.logistic.RDS')

## now the VB one
fit.sim1.vb <-
  run_vb(dat=sim1$dat, n.iter=ni, thin=th, n.chains=nc)
saveRDS(fit.sim1.vb, file='results/fit.sim1.vb.RDS')

fit.sim2.vb <-
  run_vb(dat=sim2$dat, n.iter=ni, thin=th, n.chains=nc)
saveRDS(fit.sim2.vb, file='results/fit.sim2.vb.RDS')

finish <- Sys.time()
finish-start

fit.sim1.vb <- readRDS('results/fit.sim1.vb.RDS')
fit.sim2.vb <- readRDS('results/fit.sim2.vb.RDS')
fit.sim1.logistic <- readRDS('results/fit.sim1.logistic.RDS')
fit.sim2.logistic <- readRDS('results/fit.sim2.logistic.RDS')


## Quick code to plot simulated results
sim1 <- readRDS('results/sim1.RDS')
sim2 <- readRDS('results/sim2.RDS')
post1.vb <- convert_to_df(fit.sim1.vb)
post2.vb <- convert_to_df(fit.sim2.vb)
post1.logistic <- convert_to_df(fit.sim1.logistic)
post2.logistic <- convert_to_df(fit.sim2.logistic)

## Check that abundance estimates are unbiased
png('plots/simulated_fits_abundance.png', width=8, height=6, units='in', res=500)
par(mfcol=c(2,2), mgp=c(1.5,.2, 0), tck= -.01, mar=c(3,3,1,1), oma=c(1,2,2,1))
ylim <- c(0,125)
col1 <- rgb(0,0,0,.02)
al <- c(.1,.25,.5,.75, .9)
post <- post1.logistic[,paste0("N_", 1:sim1$dat$nyrs)]
plot(0,0, type='n', xlim=c(1,ncol(post)), ylim=ylim, xlab='Year', ylab='Number Whales')
add.polygon(x=1:ncol(post), z=post, alpha.level=al)
lines(sim1$truth$abundance, col=2, lwd=2)
mtext("Logistic Model", line=1)
mtext("Data rich", side=2,line=2.5)
post <- post2.logistic[,paste0("N_", 1:sim2$dat$nyrs)]
plot(0,0, type='n', xlim=c(1,sim2$dat$nyrs), ylim=ylim, xlab='Year', ylab="Number whales")
add.polygon(x=1:ncol(post), z=post, alpha.level=al)
lines(sim2$truth$abundance, col=2, lwd=2)
mtext("Data Realistic", side=2,line=2.5)
post <- post1.vb[,paste0("N_", 1:sim1$dat$nyrs)]
plot(0,0, xlim=c(1,ncol(post)), ylim=ylim, xlab='Year', ylab="Number whales")
add.polygon(x=1:ncol(post), z=post, alpha.level=al)
lines(sim1$truth$abundance, col=2, lwd=2)
mtext("VB Model", line=1)
post <- post2.vb[,paste0("N_", 1:sim2$dat$nyrs)]
plot(0,0, xlim=c(1,sim2$dat$nyrs), ylim=ylim, xlab='Year', ylab="Number whales")
add.polygon(x=1:ncol(post), z=post, alpha.level=al)
lines(sim2$truth$abundance, col=2, lwd=2)
legend("bottom", legend=c("Truth", "Posterior"), bty='n', lty=c(1,1),
       lwd=c(2,1), col=c(2, 1))
dev.off()


## Check the VB fits for probability of availability
png('plots/simulated_fits_pa2.png', width=9, height=5, units='in', res=500)
par(mfrow=c(1,2), cex.axis=.8)
pa2.post <- post1.vb[,paste0("pa2_", 1:nseas)]
boxplot(pa2.post, names=months, ylim=c(0,1), main='Data Rich')
points(x=1:ncol(pa2.post), y=sim1$truth$pars$pa2, pch=16, col=2, cex=1.5)
pa2.post <- post2.vb[,paste0("pa2_", 1:nseas)]
boxplot(pa2.post, names=months, ylim=c(0,1), main='Data Real')
points(x=1:ncol(pa2.post), y=sim2$truth$pars$pa2, pch=16, col=2, cex=1.5)
dev.off()

png('plots/simulated_fits_pa1.png', width=9, height=5, units='in', res=500)
par(mfrow=c(1,2), cex.axis=.8)
pa1.post <- post1.vb[,paste0("pa1_", 1:sim1$dat$nyrs)]
boxplot(pa1.post, names=1:ncol(pa1.post), ylim=c(0,1), main='Data Rich')
points(x=1:ncol(pa1.post), y=sim1$truth$pars$pa1, pch=16, col=2, cex=1.5)
pa1.post <- post2.vb[,paste0("pa1_", 1:sim2$dat$nyrs)]
boxplot(pa1.post, names=1:ncol(pa1.post), ylim=c(0,1), main='Data Real')
points(x=1:ncol(pa1.post), y=sim2$truth$pars$pa1, pch=16, col=2, cex=1.5)
dev.off()

## Birth probs
png('plots/simulated_fits_zeta.png', width=9, height=5, units='in', res=500)
par(mfcol=c(2,2), cex.axis=.8, mar=c(2,.5,0,1.5), oma=c(2,5,2,.1), mgp=c(1,.5,0))
zeta.post <- post1.logistic[,paste0('zeta_', 1:sim1$dat$nyrs)]
boxplot(zeta.post, names=1:ncol(zeta.post), ylim=c(0,1))
abline(h=sim1$truth$pars$piab, col=2)
mtext("Data Rich", side=2,line=2.5)
mtext("Logistic", side=3,line=.5)
zeta.post <- post2.logistic[,paste0('zeta_', 1:sim2$dat$nyrs)]
boxplot(zeta.post, names=1:ncol(zeta.post), ylim=c(0,1))
abline(h=sim2$truth$pars$piab, col=2)
mtext("Data Real", side=2,line=2.5)
zeta.post <- post1.vb[,paste0('zeta_', 1:sim1$dat$nyrs)]
boxplot(zeta.post, names=1:ncol(zeta.post), ylim=c(0,1))
abline(h=sim1$truth$pars$piab, col=2)
mtext("VB", side=3,line=.5)
zeta.post <- post2.vb[,paste0('zeta_', 1:sim2$dat$nyrs)]
boxplot(zeta.post, names=1:ncol(zeta.post), ylim=c(0,1))
abline(h=sim2$truth$pars$piab, col=2)
mtext("Year", side=1, line=0, outer=TRUE)
dev.off()

png('plots/simulated_fits_pw.png', width=9, height=5, units='in', res=500)
par(mfrow=c(2,2), cex.axis=.8, mar=c(2.5,.5,0,1.5), oma=c(2,5,2,.1),
    mgp=c(1.5,.5,0))
bb <- 20
post <- post1.logistic[,'pw']
hist(post, xlab=NA,  main=NA, breaks=bb, xlim=c(.5,1))
mtext("Data Rich", side=3,line=.5)
mtext("Logistic", side=2,line=2.5)
post <- post2.logistic[,'pw']
hist(post, xlab=NA,  main=NA, breaks=bb, xlim=c(.5,1))
mtext("Data Real", side=3,line=.5)
post <- post1.vb[,'pw']
hist(post, xlab=NA,  main=NA, breaks=bb, xlim=c(.5,1))
mtext("VB", side=2,line=2.5)
post <- post2.vb[,'pw']
hist(post, xlab=NA,  main=NA, breaks=bb, xlim=c(.5,1))
mtext("Probability of Latent Whale (pw)", side=1, line=0, outer=TRUE)
dev.off()


png('plots/simulated_fits_avail.png', width=9, height=5, units='in', res=500)
par(mfcol=c(2,2), cex.axis=.8, mar=c(2.5,.5,0,1.5), oma=c(2,5,2,.1),
    mgp=c(1.5,.5,0))
bb <- 50
post <- post1.vb[,'mu_avail']
hist(inv.logit(post), xlab='Probability of Annual Return',  main=NA, breaks=bb, xlim=c(0,1))
mtext("Data Rich", side=2,line=2.5)
mtext("mu_avail", side=3,line=.5)
post <- post2.vb[,'mu_avail']
hist(inv.logit(post), breaks=bb, xlim=c(0,1), main=NA)
mtext("Data Real", side=2,line=2.5)
post <- post1.vb[,'sd_avail']
hist(post, main=NA, breaks=bb, xlim=c(0,1))
mtext("sd_avail", side=3,line=.5)
post <- post2.vb[,'sd_avail']
hist(post, breaks=bb, xlim=c(0,1), main=NA)
dev.off()

png('plots/simulated_fits_surv.png', width=7, height=3.5, units='in', res=500)
par(mfcol=c(1,2), cex.axis=.8, mar=c(2.5,2.5,1,1.5), oma=c(1,0,1,0),
    mgp=c(1.5,.5,0))
bb <- 20
post <- post1.vb[,'surv']
hist(post, xlab=NA,  main=NA, breaks=bb, xlim=c(.8,1))
abline(v=sim1$truth$pars$surv, col='red')
mtext("Data Rich", side=3,line=.5)
mtext("Probability of Survival", side=1,line=-.5, outer=TRUE)
post <- post2.vb[,'surv']
hist(post, xlab=NA, breaks=bb, xlim=c(.8,1), main=NA)
mtext("Data Real", side=3,line=.5)
abline(v=sim1$truth$pars$surv, col='red')
dev.off()



## ## Create matrix of pa1 by pa2 posteriors and compare that vs truth
## p.truth <- sim1$truth$pars$pa1 %o% sim1$truth$pars$pa2
## pa1.post <- post1.vb[,paste0("pa1_", 1:sim1$dat$nyrs)]
## pa2.post <- post1.vb[,paste0("pa2_", 1:sim1$dat$nseas)]
## p.post <- array(NA, dim=c(ncol(pa1.post), ncol(pa2.post), nrow(pa2.post)))
## for(i in 1:nrow(pa2.post))
##   p.post[,,i] <- as.numeric(pa1.post[i,]) %o% as.numeric(pa2.post[i,])
## matplot(p.truth-apply(p.post, 1:2, median))

## plot the distributions of effort function
for(dd in c("rich", 'real')){
  post <- if(dd=='rich') post1.vb else post2.vb
  sim <- if(dd=='rich') sim1 else sim2
  pa1.post <- post[,paste0("pa1_", 1:sim$dat$nyrs)]
  pa2.post <- post[,paste0("pa2_", 1:sim$dat$nseas)]
  k.post <- post[,'k']
  yrtemp <- 3
  ee <- seq(0,20, len=30)
  test.truth <- ldply(1:8, function(m)
    with(sim2$truth$pars,
         data.frame(month=m, effort=ee, k=k, pa1=pa1[yrtemp], pa2=pa2[m],
                    pcap=effort.fn(e=ee, pa1=pa1[yrtemp],k=k, pa2=pa2[m]))))
  test.post <- ldply(1:nrow(pa2.post), function(i)
    ldply(1:ncol(pa2.post), function(seas)
      data.frame(iter=i, effort=ee, month=seas,
                 pcap=effort.fn(ee, pa1=pa1.post[i,1],
                                pa2=pa2.post[i,seas],k=k.post[i] ))))
  g <- ggplot(test.post, aes(effort, pcap, group=iter)) + geom_line(alpha=.5) +
    facet_wrap('month') + ylim(0,1) + ylab("Probability of Detection") +
    geom_line(data=test.truth, aes(effort, pcap, group=month), lwd=.75, col='red') +
    theme_bw()
  ggsave(file=paste0('plots/simulated_fits_effort_', dd, '.png'), g, width=ggwidth, height=ggheight)
}

