## This file does Bayesian convergence checks for the model fits

### First load in the saved results and processes them where needed
message('Loading results into R..')
fit.log <- readRDS(file='results/fit.log.RDS')
fit.vb <- readRDS(file='results/fit.vb.RDS')
fit.vb_surv <- readRDS(file='results/fit.vb_surv.RDS')
fit.vb_surv_tvzeta <- readRDS(file='results/fit.vb_surv_tvzeta.RDS')
## Some convenient data.frames
post.vb <- convert_to_df(fit.vb, 'vb')
post.log <- convert_to_df(fit.log, 'log')
post.vb_surv <- convert_to_df(fit.vb_surv, 'vb_surv')
post.vb_surv_tvzeta <- convert_to_df(fit.vb_surv_tvzeta, 'vb_surv_tvzeta')

## ## Can manually check these like this
## launch_shinystan(as.shinystan(as.mcmc(fit.log)))
## launch_shinystan(as.shinystan(as.mcmc(fit.vb)))
## launch_shinystan(as.shinystan(as.mcmc(fit.vb_surv)))
## launch_shinystan(as.shinystan(as.mcmc(fit.vb_surv_tvzeta)))

### DIC comparisons
message("Calculating DIC comparisons...")
pd.vb <- fit.vb$BUGSoutput$pD
pd.log <- fit.log$BUGSoutput$pD
pd.vb_surv <- fit.vb_surv$BUGSoutput$pD
pd.vb_surv_tvzeta <- fit.vb_surv_tvzeta$BUGSoutput$pD
pds <- round(c(pd.log, pd.vb, pd.vb_surv, pd.vb_surv_tvzeta),1)
names(pds) <- c('Logistic', 'VB', 'VB Survival', 'VB Survival & Birth')
print('pD:')
print(pds)

dic.vb <- fit.vb$BUGSoutput$DIC
dic.log <- fit.log$BUGSoutput$DIC
dic.vb_surv <- fit.vb_surv$BUGSoutput$DIC
dic.vb_surv_tvzeta <- fit.vb_surv_tvzeta$BUGSoutput$DIC
dics <- c(dic.log, dic.vb, dic.vb_surv, dic.vb_surv_tvzeta)
names(dics) <- c('Logistic', 'VB', 'VB Survival', 'VB Survival & Birth')
deltaDIC <- round(dics-min(dics),1)
print('DICs:')
print(deltaDIC)

## Run through each model and store monitor values for a quick ggplot
## comparison of all models
message("Making diagnostic plots..")
diag.log <- plot_diagnostics(fit.log, 'log')
diag.vb <- plot_diagnostics(fit.vb, 'vb')
diag.vb_surv <- plot_diagnostics(fit.vb_surv, 'vb_surv')
diag.vb_surv_tvzeta <- plot_diagnostics(fit.vb_surv_tvzeta, 'vb_surv_tvzeta')
diag.all <- rbind(diag.log, diag.vb, diag.vb_surv, diag.vb_surv_tvzeta)
diag.all.long <- melt(diag.all, id.vars=c("model", "par"))
g <- ggplot(diag.all.long, aes(model, value, color=model)) +
  geom_jitter(width=.15, height=0, alpha=.5) +
  facet_wrap("variable", scales='free_y', ncol=1) + theme_bw()
ggsave('plots/fit_diag_all.png', g, width=ggwidth, height=ggheight)

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

message("Plotting posterior predictive distributions...")
pp.log <- get_post_predict(post.log)
pp.vb <- get_post_predict(post.vb)
pp.vb_surv <- get_post_predict(post.vb_surv)
pp.vb_surv_tvzeta <- get_post_predict(post.vb_surv_tvzeta)
## Combine together
pp <- rbind(cbind(model='vb', pp.vb),
            cbind(model='vb_surv', pp.vb_surv),
            cbind(model='vb_surv_tvzeta', pp.vb_surv_tvzeta),
            cbind(model='log', pp.log))

## plot all on the same figure against the data
yy <- data.frame(year=years, whales_observed=apply(apply(dat$x,1:2, function(p) sum(p)>0), 2, sum))
g <- ggplot(pp, aes(year, whales_observed)) +
  geom_jitter(height=.1, width=.5, size=.1, alpha=.1) +
  facet_wrap('model')
g <- g+ geom_point(data=yy, col='red', size=3)
ggsave('plots/posterior_predictive.png', g, width=ggwidth, height=ggheight)

## Compare the posterior predictive with and without the tvzeta component.
pp2 <- droplevels(subset(pp, model %in% c('vb_surv','vb_surv_tvzeta')))
ggplot(pp2, aes(x=factor(year), y=whales_observed, fill=model))+
  geom_violin(scale='width')
