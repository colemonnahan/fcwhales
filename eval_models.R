## This file does Bayesian convergence checks for the model fits

## DIC comparisons
dic.vb1 <- fit.vb1$BUGSoutput$DIC
dic.log1 <- fit.log1$BUGSoutput$DIC
dic.vb2 <- fit.vb2$BUGSoutput$DIC
dic.log2 <- fit.log2$BUGSoutput$DIC
dic.vb_surv <- fit.vb_surv$BUGSoutput$DIC
dics <- rbind(c(dic.vb1, dic.log1), c(dic.vb2, dic.log2))
colnames(dics) <- c("VB", "logistic")
row.names(dics) <- c("Pared Data", "Full Data")
print("DIC table")
print(dics)
print(dic.vb_surv)

## Function to automatically make some diagnostic plots for a given fitted
## model.
plot_diagnostics <- function(fit, name){
  N <- nrow(fit$BUGSoutput$sims.list[[1]])
  fit <- as.mcmc(fit)
  pdf(file=paste0("plots/fit_diag_", name, ".pdf"), width=8, height=6,
      onefile=TRUE)
  mon <- data.frame(monitor(convert_array(fit), warmup=0, probs=.5, print=FALSE))
  mon$par <- rownames(mon)
  mon$pct.ess <- 100*mon$n_eff/N
  par(mfrow=c(2,1))
  x1 <- mon[order(mon$pct.ess)[1:10], ]
  barplot(x1$pct.ess, names=x1$par, main='% ESS', ylim=c(0,100));box()
  x2 <- mon[order(mon$Rhat, decreasing=TRUE)[1:10], c('par', 'Rhat')]
  barplot(x2$Rhat-1, names=x2$par, main='Rhat - 1.0', ylim=c(0,.5));box()
  abline(h=.1)
  par(mfrow=c(3,3))
  ## geweke.plot(fit)
  traceplot(fit)
  dev.off()
  mon$model <- name
  mon <- subset(mon, select=c("model", "pct.ess", 'n_eff', "Rhat", "par"))
  row.names(mon) <- NULL
  return(mon)
}

## Run through each model and store monitor values for a quick ggplot
## comparison of all models
diag.vb1 <- plot_diagnostics(fit.vb1, 'vb1')
diag.vb2 <- plot_diagnostics(fit.vb2, 'vb2')
diag.log1 <- plot_diagnostics(fit.log1, 'log1')
diag.log2 <- plot_diagnostics(fit.log2, 'log2')
diag.vb_surv <- plot_diagnostics(fit.vb_surv, 'vb_surv')
diag.all <- rbind(diag.vb1, diag.vb2, diag.log1, diag.log2, diag.vb_surv)
diag.all.long <- melt(diag.all, id.vars=c("model", "par"))
g <- ggplot(diag.all.long, aes(model, value, color=model)) +
  geom_jitter(width=.15, height=0, alpha=.5) +
  facet_wrap("variable", scales='free') + theme_bw()
ggsave('plots/fit_diag_all.png', g, width=ggwidth, height=ggheight)

 ## launch_shinystan(as.shinystan(fit.vb1))
 ## launch_shinystan(as.shinystan(fit.vb2))
## launch_shinystan(as.shinystan(fit.log1))
## launch_shinystan(as.shinystan(fit.log2))

### As suggested by a reviewer we compare posterior predictive
### distributions between model versions. That is, inside the model we
### simulate data given the number of whales alive and the probability of
### observing them (one for each iteration of MCMC to get a distribution)
### and then compare that to the real data. The real data should reasonably
### come from these distributions. See Gelman et al. (2011) section 6.3).

get_post_predict <- function(x){
### Quick function to extract the results and plot them with ggplot, where
### 'x' is the posteriors as data.frames from load_results
  xx <- x[, paste('post_predict', 1:nyrs, sep='_')]
  names(xx) <- 2003+1:nyrs
  temp <- melt(xx, variable.name='year', value.name='whales_observed', id.vars=NULL)
  temp$year <- as.numeric(as.character(temp$year))
  return(temp)
}

pp.vb1 <- get_post_predict(post.vb1)
pp.vb_surv <- get_post_predict(post.vb_surv)
pp.log1 <- get_post_predict(post.log1)
## Combine together
pp <- rbind(cbind(model='vb1', pp.vb1), cbind(model='vb_surv', pp.vb_surv),
            cbind(model='log1', pp.log1))

## plot all on the same figure
yy <- data.frame(year=2003+1:nyrs, whales_observed=apply(apply(dat2$x,1:2, function(p) sum(p)>0), 2, sum))
g <- ggplot(pp, aes(year, whales_observed)) + geom_jitter() + facet_wrap('model')
print(g+ geom_point(data=yy, col='red', size=3))
