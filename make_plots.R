library(reshape2)
library(plyr)

### Explore data ---------------
## Scatterplot of effort vs unique whales. Sort of approximates
## pcap=f(effort). Note I set ylim to (0,85) since there's about 85 whales
## from the model. So gives sense of the proportion captured
g <- ggplot(obs.effort, aes(effort, unique.whales)) +
  geom_point(aes( color=factor(year), shape=year>2015), size=2, alpha=.5) +
  ylim(0,85) + facet_wrap('month') + ylab("Unique whales observed") + theme_bw()
ggsave('plots/captures_by_effort.png', g, width=ggwidth,
       height=ggheight)

whales.table.long <- melt(whales.table, id.vars='year')
write.table(whales.table, file='whales.table.csv', sep=',', row.names=FALSE)
g <- ggplot(whales.table.long, aes(year, value, group=variable, color=variable)) + geom_line()
ggsave('plots/captures_by_year.png', g, width=ggwidth,
       height=ggheight)

unique.whales.table <- dcast(ddply(obs.effort, .(year, month), summarize,
            out=paste0(unique.whales, " (", effort, ")")), value.var='out', year~month)
write.table(unique.whales.table, file='unique.whales.table.csv', sep=',', row.names=FALSE)

g <- ggplot(obs.effort, aes(month, unique.whales/effort)) +
  geom_jitter(aes( color=year %in% c(2007, 2012)), size=2, height=0,
  width=.1, alpha=.6) + ylim(0,9) + ylab("(Unique whales observed) / (Days Effort)")
ggsave('plots/relative_captures_by_month.png', g, width=ggwidth,
       height=ggheight)

## Effort vs unique whales observed
obs.effort.long <- melt(obs.effort, id.vars=c('month', 'year'))
g <- ggplot(obs.effort.long, aes(month, value, group=variable, color=variable)) + geom_line() +
  facet_wrap("year", scales='fixed')
ggsave('plots/effort_observations_by_year.png', g, width=10, height=ggheight)

## Explore individual whale seasonality. Need to be careful with effort
## since it is unbalanced. Idea is to check whether some whales are only
## ever seen in certain months. This would be a problem for the model but
## more importantly for designing surveys.
x.long <- melt(x)
## Which whale IDs were seen in > 3 years
temp <- ddply(x.long, .(ID, year), summarize, months.obs=sum(value==1))
## Note some of these can be 0 since we took out some years from the data,
## and there are augmented "unseen" whales.
annual.counts <- ddply(temp, .(ID), summarize, years.obs=sum(months.obs>0))
annual.counts <- subset(annual.counts, years.obs>0)
IDs <- subset(annual.counts, years.obs > 7)$ID
xx <- merge(x.long, obs.effort)
xx <- subset(xx, ID %in% IDs & effort > 0)
## This shows the annual (rough) capture probabilities for individual
## whales which were seen more than 7 years. The idea is to make sure
## whales aren't only there in certain months.
g <- ggplot(xx, aes(month, value/effort)) + geom_point(aes(size=effort), alpha=.25) +
  geom_line(aes(group=year), alpha=.25) + facet_wrap('ID')  +
  scale_size(range = c(0, 3)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90)) + xlab("Month") +
  ylab("Observed/Days Effort")
ggsave('plots/individual_sightings.png', g, width=10, height=1.5*ggheight)

## Which whales were seen before last year
IDs <- subset(ddply(subset(temp, months.obs>0), .(ID), summarize,
             year.obs=min(year)), year.obs < 2016)$ID
## which whales were seen in  > 1 year, tells us about resighting
temp2 <- ddply(subset(x.long, ID %in% IDs), .(ID, year), summarize, times=sum(value))
temp3 <- ddply(temp2, .(ID), summarize, times=sum(times>0))
##table(temp3$times)
print("Percentage of whales seen before 2016 that were seen > 1 time")
print(round(100*mean(temp3$times > 1),2))


### Explore model fits ---------------
## Compare some key estimates between model combinations
post <- rbind.fill(post.vb1, post.vb2, post.log1, post.log2)
post$model <- as.factor(mapvalues(post$model, from=c('log1', 'log2','vb1', 'vb2'),
                        to=c('Logistic', 'Logistic2', 'VonBert',
                     'VonBert2')) )
post.vb <- droplevels(subset(post, model %in% c('VonBert', 'VonBert2')))
post.log <- droplevels(subset(post, ! model %in% c('VonBert', 'VonBert2')))

g <- ggplot(post) + geom_violin(aes(model, surv)) + ylab("Survival") + ylim(.8,1)
ggsave('plots/posterior_surv.png', g, width=ggwidth, height=ggheight)
g <- ggplot(post) + geom_violin(aes(model, TotalN)) + ylab("Total Whales")
ggsave('plots/posterior_totalN.png', g, width=ggwidth, height=ggheight)

## Effort relationship
g <- ggplot(post.vb) + geom_violin(aes(model, k))
ggsave('plots/posterior_k.png', g, width=ggwidth, height=ggheight)

## Show priors for availability hyperparameters
xx <- post.vb[, c('model', 'mu_avail', 'sd_avail')]
xx <- melt(xx, id.vars='model')
## mu_avail ~ dnorm(2,1/(.1*.1))
## sd_avail ~ dnorm(0,10)T(0,) # half-normal
xseq <- seq(-3,3, len=100)
p1 <- data.frame(variable='mu_avail', x=xseq, y=dnorm(xseq, 0,1))
xseq <- seq(0, 3, len=100)
p2<- data.frame(variable='sd_avail', x=xseq, y=dnorm(xseq, 0, 2))
p3 <- data.frame(model=c('VonBert', "VonBert2"),rbind(p1,p2))
g <- ggplot(xx, aes(x=value)) + geom_histogram(aes(y=..density..),
                                          position="identity") +
  geom_line(data=p3, aes(x,y))+
  facet_grid(model~variable, scales='free')
ggsave('plots/posterior_avail.png', g, width=ggwidth, height=ggheight)


xx <- post.vb[, c('model',paste0('pa2_', 1:dat$nseas))]
xx <- melt(xx, id.vars='model')
xx$month <- mapvalues(xx$variable, from=paste0('pa2_', 1:dat$nseas), to=months)
pa2.meds <- ddply(subset(xx, model='VonBert'), "month", summarise, median = median(value))
g <- ggplot(xx)+ geom_violin(aes(x=month, value)) + facet_wrap("model", nrow=2) +
  ylab("Monthly Availability") + ylim(0,1)
ggsave('plots/posterior_pa2.png', g, width=ggwidth, height=ggheight)

xx <- post.vb[, c('model',paste0('pa1_', 1:dat$nyrs))]
xx <- melt(xx, id.vars='model')
xx <- na.omit(xx)
xx$year <- -1 + ifelse(xx$model=='VonBert2', 2002, 2004) + as.numeric(substr(as.character(xx$variable),5,7))
pa1.meds <-
  ddply(subset(xx, model='VonBert'), "year", summarise, median = median(value))
g <- ggplot(xx)+ geom_violin(aes(x=factor(year), value)) +
  facet_wrap("model", nrow=2) + ylab("Annual Availability")
ggsave('plots/posterior_pa1.png', g, width=ggwidth, height=ggheight)

xx <- post.log[, c('model','gamma2',paste0('gamma3_', 1:dat$nseas))]
xx <- melt(xx, id.vars=c('model', 'gamma2'))
xx$month <- mapvalues(xx$variable, from=paste0('gamma3_', 1:dat$nseas), to=months)
xx$pcap <- inv.logit(xx$value+xx$gamma2*10)
g <- ggplot(xx)+ geom_violin(aes(x=month, pcap)) + facet_wrap("model", nrow=2)
ggsave('plots/posterior_pcap.png', g, width=ggwidth, height=ggheight)


## ## Effort vs data
## pa1.meds <-
##   ddply(subset(xx, model='VonBert'), "year", summarise, median = median(value))
## ggplot(obs.effort, aes(effort, unique.whales)) +
##   geom_point(aes( color=factor(year), shape=year<2004), size=2) + ylim(0,85) + facet_wrap('month') + geom_smooth()

## Population trends
obs.df <- ldply(c("VonBert", "VonBert2", "Logistic", "Logistic2"),
                function(x) data.frame(model=x, year=years, y=whales.table$unique.whales))
xx <- post[, c('model',paste0('N_', 1:dat$nyrs))]
xx <- na.omit(melt(xx, id.vars='model'))
xx$year <- -1 + ifelse(xx$model %in% c('VonBert2', 'Logistic2'), 2002, 2004) + as.numeric(substr(as.character(xx$variable),3,5))
xx <- ddply(xx, .(model, year), summarize,
            lwr=quantile(value, .025),
            med=median(value),
            upr=quantile(value, .975))
g <- ggplot(xx)+ geom_ribbon(aes(x=year, ymin=lwr, ymax=upr)) +
  facet_wrap("model") + geom_line(aes(x=year, y=med), col='red', lwd=1) +
  ylim(0,110) + ylab("Abundance") + geom_line(data=obs.df, aes(year, y=y), lty=2)
ggsave('plots/posterior_abundance.png', g, width=ggwidth, height=ggheight)

### End of final code
stop("dont source past this point")


## #### The rest is temporary or old or unused and probably wont work, but is
## #### left here since it may be useful

## ### Explore priors vs posterior ---------------
## ## Rough estimates of pcap
## ind.vec <- grep(x=varnames(mod1),'gamma3')
## post <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])

## ### NEED TO CHECK JAGS dt function below
## col.post <- gray(.5)
## lty.prior <- 2; col.prior <- 1
## png('plots/priors_vs_posteriors1.png', units='in', width=3,height=6, res=500)
## par(mfrow = c(3,1))
## #priors for beta dt(0,0.16,3)
## ind <- which(varnames(mod1) =="beta")
## post <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
## an.x <- seq(.9*min(post), 1.1*max(post), len=1000)
## hist(x=post, xlab=NA , main = "beta", ylab = "density", freq=FALSE, col=col.post)
## lines(an.x, dt(an.x, df=3/0.16), lty=lty.prior, col=col.prior)
## ## priors for gamma0 dnorm(0, 0.3)
## ind <- which(varnames(mod1) =="gamma2")
## post <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
## an.x <- seq(.9*min(post), 1.1*max(post), len=1000)
## hist(x=post, xlab=NA , main = "gamma2", ylab = "density", freq=FALSE, col=col.post)
## lines(an.x, dnorm(an.x, 0, sd = 1.82), lty=lty.prior, col=col.prior)
## #priors for sd ~ U(0,10)
## ind <- which(varnames(mod1) =="sd")
## post <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
## an.x <- seq(.9*min(post), 1.1*max(post), len=1000)
## hist(x=post, xlab=NA , main = "sd", ylab = "density", freq=FALSE, col=col.post)
## lines(an.x, dunif(an.x,0,10), lty=lty.prior, col=col.prior)
## dev.off()
## ## priors for gamma0 dnorm(0, 0.3)
## png('plots/priors_vs_posteriors2.png', units='in', width=6,height=6, res=500)
## par(mfrow=c(3,4))
## for(i in 1:(nyrs-1)){
## par <- paste0("zeta[", i,"]")
## ind <- which(varnames(mod1) ==par)
## post <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
## an.x <- seq(0,1, len=1000)
## hist(x=post, xlab=NA , xlim=c(0,1), main = par, ylab = "density", freq=FALSE, col=col.post)
## lines(an.x, dbeta(an.x, 1, nyrs-i), lty=lty.prior, col=col.prior)
## }
## dev.off()
## ## priors for gamma3 dt(0,0.16, 3)
## png('plots/priors_vs_posteriors3.png', units='in', width=6,height=6, res=500)
## par(mfrow=c(3,3))
## for(i in 1:nseas){
## par <- paste0("gamma3[", i,"]")
## ind <- which(varnames(mod1) ==par)
## post <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
## an.x <- seq(0,1, len=1000)
## hist(x=post, xlab=NA ,  main = par, ylab = "density", freq=FALSE, col=col.post)
## lines(an.x, dt(an.x, df=3, ncp=.16), lty=lty.prior, col=col.prior)
## }
## dev.off()

## ## Posterior for pcap by year/seas
## ## model=
## ### logit(p[ind,yr,seas]) <-
## ### gamma3[seas] + gamma2*(effort[yr,seas]-effort.bar)

## effort.bar <- mean(effort)
## pcap.fn.seas <- NULL
## pcap.post <- array(NA, dim=c(nyrs, nseas, 3),
##                     dimnames=list("year"=years, "month"=months,
##                                   "quantile"=c("lwr", "med", "upr")))
## for(yr in 1:nyrs){
##   ind <- which(varnames(mod1) == 'gamma2')
##   gamma2 <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
##   for(seas in 1:nseas){
##     par <- paste0("gamma3[", seas,"]")
##     ind <- which(varnames(mod1) ==par)
##     gamma3 <- c(mod1[[3]][,ind], mod1[[1]][,ind], mod1[[2]][,ind])
##     ## explore shape as function of effort
##     effort.vec <- seq(-10, 20, len=10)
##     if(yr==1) ## technically shoudl do this for all years
##       pcap.fn.seas[[seas]] <- ldply(1:length(effort.vec), function(i)
##         data.frame(season=seas, effort=effort.vec[i],
##                    pcap=inv.logit(x=gamma3 + gamma2*(effort.vec[i]-effort.bar))))
##     ## The one used in making predictions in the model
##     pcap.post[yr,seas,] <-quantile(x=inv.logit(x=gamma3 +
##                                                  gamma2*(effort[yr,seas]-effort.bar)),
##                                    probs=c(0.025, .5, .975))
##   }
## }

## pcap.long <- melt(pcap.post)
## pcap.long <- dcast(pcap.long, year+month~quantile, value.var='value')
## g <- ggplot(pcap.long, aes(year, month, color=upr, size=med)) + geom_point()
## g

## pcap2 <- ddply(do.call(rbind, pcap.fn.seas), .(effort, season), summarize,
##                lwr=quantile(pcap, 0.025),
##                med=quantile(pcap, 0.5),
##                upr=quantile(pcap, 0.975))
## g <- ggplot(pcap2, aes(effort, med, fill=factor(season))) + ylim(0,1)+
##                geom_ribbon(aes(ymin=lwr, ymax=upr)) + facet_wrap("season")

## pcap3 <- subset(do.call(rbind, pcap.fn.seas), effort==0)
## g <- ggplot(pcap3, aes(factor(season), pcap)) + geom_violin()
## g


## ## Quickly explore some effort shapes
## ff <- function(gammat, gamma1, eff) 1/(1+exp(-gammat - gamma1*eff))

## res <- ldply(c(1,2), function(gamma1)
##   ldply(c(-20, -10), function(gammat)
##     ldply(seq(0,15, len=50), function(eff)
##       data.frame(gamma1=gamma1, gammat=gammat, effort=eff, pcap=ff(gammat,
##                                                                    gamma1, eff))
##       )))
## ggplot(res, aes(effort, pcap, group=gamma1)) + geom_line() + facet_wrap("gammat")


## ef <- function(effort, gamma, pseas, pyear)
##   pseas*pyear*(1-exp(-gamma*effort))
## pyears <- runif(5, .7, 1)
## pseas <- c(seq(.01, .4, len=4), seq(.45, .05, len=4))
## res <- ldply(1:5, function(yr)
##   ldply(1:8, function(seas){
##     ee <- seq(0, 20, len=50)
##     data.frame(year=yr, month=seas, effort=ee, pcap=ef(ee, .2, pseas[seas], pyears[yr]))
##  }))
## ggplot(res, aes(effort, pcap, group=month, color=factor(month))) + geom_line() + facet_wrap("year")

## ## Explore prior for random effects for probabilities

## ## hyperdistributions and randomly generated values, trying to get it close
## ## to uniform... from this is looks like mu~N(0,1) and sd~N(0,2)T(0,) gives
## ## uniform random effects a priori which seems good.
## mu <- rnorm(5000, 0, 1)
## sd.mu <- abs(rnorm(5000, 0,2))
## pp <- data.frame(mu=mu, sd=sd.mu, prob=inv.logit(rnorm(n=5000, mean=mu, sd.mu)))
## ggplot(pp, aes(prob)) + geom_histogram() + xlim(0,1)
## ## ggplot(pp, aes(logit(prob))) + geom_histogram()

## ## Explore hyper-prior on annual random effects
## mu <- rnorm(5000, 2.2, .05)
## sd.mu <- abs(runif(5000, 0,7))
## pp <- data.frame(mu=mu, sd=sd.mu, prob=inv.logit(rnorm(n=5000, mean=mu, sd.mu)))
## ggplot(pp, aes(prob)) + geom_histogram() + xlim(0,1)


## ### Explore model fits

## #Evaluate model estimates:
## #list of variable names in the model
## varnames(mod1)

## #Parameter estimates:
## gamma2<-summary(mod1[,which(varnames(mod1) =="gamma2")], quantiles = c(0.5, 0.025, 0.975))[[2]]
## #gamma1<-summary(mod1[,which(varnames(mod1) =="gamma1")], quantiles = c(0.5, 0.025, 0.975))[[2]]

## #monthly rates of recapture
## gamma3<- summary(mod1[,which(varnames(mod1) =="gamma3[1]"):which(varnames(mod1) == "gamma3[8]"), ], quantiles = c(0.5, 0.025, 0.975) )[[2]]

## surv.est<-inv.logit(summary(mod1[,which(varnames(mod1) =="beta")], quantiles = c(0.5, 0.025, 0.975) )[[2]] )

## #Abundance estimates:
## est.N<- summary(mod1[,which(varnames(mod1) =="N[1]"):which(varnames(mod1) == "N[13]"), ], quantiles = c(0.5, 0.025, 0.975) )[[2]]

## #true.n<-dget("True.Abundance.txt")
## #true.n<-n.whales
## #FIRST PLOT
## plot(2002:2014, est.N[,1], type = 'l', xlab = "Time", ylab = "Abundance", ylim = c(0,120))
## lines(2002:2014, est.N[,2], lty = 2)
## lines(2002:2014, est.N[,3], lty = 2)

## #SECOND PLOT - plot the summer years
## years<-2003:2015
## plot(years, est.N[,1], pch=15, xlab = "Year (Spring)", ylab = "Abundance", ylim = c(0,120))
## for(i in 1:k){
##   lines(c(years[i], years[i]), c(est.N[i,2], est.N[i,3]) )
## }

## #Primera figura - el verano
## years<-2003:2015
## plot(years, est.N[,1], pch=15, xlab = "Año (verano)", ylab = "Abundancia", ylim = c(0,120))
## for(i in 1:k){
##   lines(c(years[i], years[i]), c(est.N[i,2], est.N[i,3]) )
## }

## #en ingles:
## #Primera figura - el verano
## years<-2003:2012
## plot(years, est.N[,1], pch=15, xlab = "Year (summer)", ylab = "Abundance", ylim = c(0,100))
## for(i in 1:k[1]){
##   lines(c(years[i], years[i]), c(est.N[i,2], est.N[i,3]) )
## }


## write.table(est.N, "abundancia.csv", sep = ',')

## ###  ESTIMATE THE GROWTH RATE FOR THE 2003-2015 PERIOD ########

## #obtain matrix of abundance estimates from model:
## N1.ind<-which(varnames(mod1) == "N[1]")
## N13.ind<- which(varnames(mod1) == "N[13]")
## abund.vals<- append( mod1[[3]][,N1.ind:N13.ind], append(mod1[[1]][, N1.ind:N13.ind], mod1[[2]][,ind]) )

## abund.vals<- rbind( mod1[[3]][,N1.ind:N13.ind], mod1[[1]][, N1.ind:N13.ind])

## years1<- 2003:2012
## an.x1<-rep(1, times = length(years1))
## an.x2<-1:length(years1)
## an.x<-rbind(an.x1, an.x2)
## an.xtx<-an.x%*%t(an.x)
## inv.an.x<-solve(an.xtx)

## b.end<-matrix(ncol=2,nrow=length(abund.vals[,1]) )
## for( i in 1:length(abund.vals[,1])){
## 	b.end[i,]<-t( inv.an.x%*%an.x%*%(log(abund.vals[i,1:10]) ) )
## 	}
## b.end.q<-quantile(b.end[,2], probs = c(0.5, 0.025, 0.975) )


## #TREND FROM 2004 TO 2012
## years2<- 2005:2015
## an.x1<-rep(1, times = length(years2))
## an.x2<-1:length(years2)
## an.x<-rbind(an.x1, an.x2)
## an.xtx<-an.x%*%t(an.x)
## inv.an.x<-solve(an.xtx)

## b.end2<-matrix(ncol=2,nrow=length(abund.vals[,1]) )
## for( i in 1:length(abund.vals[,1])){
## 	b.end2[i,]<-t( inv.an.x%*%an.x%*%(log(abund.vals[i,3:13]) ) )
## 	}
## b.end.q2<-quantile(b.end2[,2], probs = c(0.5, 0.025, 0.975) )


## years3<- 2005:2014
## an.x1<-rep(1, times = length(years3))
## an.x2<-1:length(years3)
## an.x<-rbind(an.x1, an.x2)
## an.xtx<-an.x%*%t(an.x)
## inv.an.x<-solve(an.xtx)

## b.end3<-matrix(ncol=2,nrow=length(abund.vals[,1]) )
## for( i in 1:length(abund.vals[,1])){
## 	b.end3[i,]<-t( inv.an.x%*%an.x%*%(log(abund.vals[i,3:12]) ) )
## 	}
## b.end.q3<-quantile(b.end3[,2], probs = c(0.5, 0.025, 0.975) )


## #PROBABILITY OF CAPTURE RELATIONSHIPS:

## effort.bar<- mean(effort)
## new.effort<- 1:11
## gamma2.ind<- which(varnames(mod1) =="gamma2")
## gam.eff<- append( mod1[[3]][,gamma2.ind], append(mod1[[1]][,gamma2.ind], mod1[[2]][,gamma2.ind]) )

## gam3.jan.ind<- which(varnames(mod1) =="gamma3[3]")
## gam.jan<- append( mod1[[3]][,gam3.jan.ind], append(mod1[[1]][,gam3.jan.ind], mod1[[2]][,gam3.jan.ind]) )

## pr.new<- matrix(ncol = 11, nrow = 6000)
## pr.new.q<- matrix(ncol = 11, nrow = 3)

## inv.logit(mean(gam.jan) + mean(gam.eff)*(new.effort - effort.bar) )

## for(i in 1:length(new.effort)){
## pr.new[,i]<- inv.logit(gam.jan +  gam.eff*(new.effort[i] - effort.bar)  )
## pr.new.q[,i]<- quantile(pr.new[,i], probs = c(0.5, 0.025, 0.975) )
## }

## par(mfrow = c(1,2))
## #  MONTH PLOT
## plot(1:8, inv.logit(gamma3[,1]) , pch = 15, xlab = "Month",  ylab = "Pr(capture)" , xaxt = 'n', ylim = c(0, 0.32))
## axis(1, at = 1:8, labels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"), las = 2)
## for(i in 1:8){
## 	lines(c(i, i), c(inv.logit(gamma3[i,2]), inv.logit(gamma3[i,3])  ) )
## }
## text(1.5, 0.3, "A")
## #EFFORT PLOT
## plot(new.effort, pr.new.q[1,] , pch = 15, xlab = "Effort(days)",  ylab = "Pr(capture)" , ylim = c(0.1, 0.65) )
## for(i in 1:11){
## 	lines(c(i, i), c(pr.new.q[2,i], pr.new.q[3,i])  )
## }
## text(1.5, 0.61, "B")


## #PLOT OF POSTERIOR DISTRIBUTIONS WITH PRIORS
## varnames(mod1[,10:15])
## library(KernSmooth)
## par(mfrow = c(2,2))

## #priors for beta dt(0,0.16,3)
## an.x<-seq(-15, 15, by = 0.1)
## plot(an.x, dt(an.x, df=(3/0.16), ncp=0), type = 'l', lty=2, xlab = "beta", ylab = "density" , ylim = c(0, 3.5))
## beta.ind<- which(varnames(mod1) =="beta")
## lines( bkde( append( mod1[[3]][,beta.ind], append(mod1[[1]][,beta.ind], mod1[[2]][,beta.ind]) )  )  )

## #priors for gamma0 dnorm(0, 0.3)
## an.x<-seq(-15, 15, by = 0.1)
## plot(an.x, dnorm(an.x, 0, sd = 1.82), type = 'l', lty=2, xlab = "gamma", ylab = "density", ylim = c(0,0.82) )
## gamma.ind<- which(varnames(mod1) =="gamma2")
## lines( bkde( append( mod1[[3]][,gamma.ind], append(mod1[[1]][,gamma.ind], mod1[[2]][,gamma.ind]) )  )  )
## post <- c(mod1[[3]][,gamma.ind], mod1[[1]][,gamma.ind], mod1[[2]][,gamma.ind])
## hist(x=post , xlab = "gamma", ylab = "density", freq=FALSE)
## lines(an.x, dnorm(an.x, 0, sd = 1.82), lty=2)

## ##priors for gamma1 dnorm(0, 0.3)
## #an.x<-seq(-15, 15, by = 0.1)
## #plot(an.x, dnorm(an.x, 0, sd = 1.82), type = 'l', lty=2, xlab = "gamma 1", ylab = "density", ylim = c(0,0.65) )
## #lines( bkde( append( mod1[[3]][,12], append(mod1[[1]][,12], mod1[[2]][,12]) )  )  )

## #priors for sd ~ U(0,10)
## an.x<-seq(0, 10, by = 0.5)
## plot(an.x, dunif(an.x, 0,10), type = 'l', lty=2, xlab = "sd", ylab = "density", ylim = c(0,0.8) )
## sd2.ind<- which(varnames(mod1) =="sd")
## lines( bkde( append( mod1[[3]][,sd2.ind], append(mod1[[1]][,sd2.ind], mod1[[2]][,sd2.ind]) )  )  )

## #no random effects on the pr. of survival for this model
## if(F){
## #priors for sd ~ U(0,2)
## an.x<-seq(0, 10, by = 0.5)
## plot(an.x, dunif(an.x, 0,2), type = 'l', lty=2, xlab = "sd 1", ylab = "density", ylim = c(0,0.8) )
## sd1.ind<- which(varnames(mod1) =="sd[1]")
## lines( bkde( append( mod1[[3]][,sd1.ind], append(mod1[[1]][,sd1.ind], mod1[[2]][,sd1.ind]) )  )  )
## }


## #CROSS CORRELATION OF MODEL PARAMETERS
## par(mfrow = c(1,1))
## crosscorr(mod1[,c(10:12,14:15)])
## crosscorr.plot(mod1[,c(10:12,14:15)], col = gray(level = seq(0,1, by = 0.20) ), cex.lab = 0.75, cex.axis = 0.75 )
## crosscorr.plot(mod1, col = gray(level = seq(0,1, by = 0.20) ), cex.lab = 0.75, cex.axis = 0.75, cex = 0.5)

## par(mfrow=c(1,1))
## #Total Number of whales
## last.N.ind<-which(varnames(mod1) == "N[13]")
## last.N<- append( mod1[[3]][,last.N.ind], append(mod1[[1]][,last.N.ind], mod1[[2]][,last.N.ind]) )

## lambda.ind<-which(varnames(mod1) == "lambda")
## the.lambda<- append( mod1[[3]][,lambda.ind], append(mod1[[1]][,lambda.ind], mod1[[2]][,lambda.ind]) )

## par(mfrow = c(1,2))
## hist(the.lambda, xlab = "All Whales 2003 - 2015", prob = T, breaks = 160:200, main = '', ylim = c(0,0.12), col = "gray", lwd = 1)
## hist(last.N , add = F, breaks = 80:110, prob=T, lwd = 1, xlab = "Abundance 2015", main = '')
## #legend(20, 0.3, c("Abundance 2015", "All Whales 2003-2015"), col = c("black", "gray"), lty=c(2,1), lwd = c(3,3) )


## ## Figura 3 -
## hist(the.lambda, xlab = "Abundancia", prob = T, breaks = seq(70,190, by = 2) , main = '', ylim = c(0,0.1), col = "gray", lwd = 3, ylab = "Densidad de probabilidad")
## hist(last.N , add = T, breaks = seq(70,190, by = 2) , prob=T, lty = 2, lwd = 3)
## legend(80, 0.1, c("Abundancia 2012", "Total 2003 - 2012"), col = c("black", "gray"), lty=c(2,1), lwd = c(3,3) )

## last.N.q<- quantile(last.N, probs = c(0.5, 0.025, 0.975) )
## the.lambda.q<-quantile(the.lambda, probs = c(0.5, 0.025, 0.975) )

## #ANOTHER PLOT - NOT USED
## plot( bkde( last.N) , type = 'l', xlab = "Whale Abundance", ylab = "density", ylim = c(0,0.26), xlim=c(15,60), lty = 2  )

## lines(bkde( append( mod1[[3]][,9], append(mod1[[1]][,9], mod1[[2]][,9]) )  ) )

## #SURVIVAL RATE BY YEAR - no  etas
## beta.ind<-which(varnames(mod1) == "beta")
## the.beta<-append( mod1[[3]][,beta.ind], append(mod1[[1]][,beta.ind], mod1[[2]][,beta.ind]) )
## #the.etas<-matrix(nrow = length(the.beta), ncol = k)
## #eta.ind<-which(varnames(mod1) == "etas[1]")
## #for(i in 1:k){
## #	the.etas[,i]<-append( mod1[[3]][,eta.ind + (i-1)], append(mod1[[1]][,eta.ind + (i-1)], mod1[[2]][,eta.ind + (i-1)]) )
## #}
## an.surv<- inv.logit(the.beta)# + the.etas)
## an.surv.q<-matrix(nrow=3,ncol=1)

## an.surv.q<-quantile(an.surv, probs = c(0.5, 0.025, 0.975) )

## plot(bkde(an.surv), type = 'l', xlab = "Probability of survival", ylab = "Density")

## #not using random effects on survival
## if(F){
## plot(years, an.surv.q[1,], pch = 15, xlab = "Year (Spring)", ylab = "Survival", ylim = c(0.5, 1))
## for(i in 1:k){
## 	lines(c(years[i], years[i]), c(an.surv.q[2,i], an.surv.q[3,i]) )
## }

## library(vioplot)
## vioplot(an.surv[,1], an.surv[,2], an.surv[,3], an.surv[,4], an.surv[,5], an.surv[,6], an.surv[,7], an.surv[,8], an.surv[,9],names = c("2002", "2003", "2004", "2005","2006", "2007", "2008", "2009", "2010" ), col = "gray")
## }

## #ANNUAL CAPTURE RATE BY YEAR
## gamma.ind<-which(varnames(mod1) == "gamma")
## the.gamma<-append( mod1[[3]][,gamma.ind], append(mod1[[1]][,gamma.ind], mod1[[2]][,gamma.ind]) )

## the.etap<-matrix(nrow = length(the.gamma), ncol = k)
## etap.ind<-which(varnames(mod1) == "etap[1]")
## for(i in 1:k){
## 	the.etap[,i]<-append( mod1[[3]][,etap.ind + (i-1)], append(mod1[[1]][,etap.ind + (i-1)], mod1[[2]][,etap.ind + (i-1)]) )
## }
## an.capt<- inv.logit(the.gamma + the.etap)
## an.capt.q<-matrix(nrow=3,ncol=k)
## for(i in 1:k){
## 	an.capt.q[,i]<-quantile(an.capt[,i], probs = c(0.5, 0.025, 0.975) )
## }

## plot(years, an.capt.q[1,], pch = 15, xlab = "Year (Summer)", ylab = "Pr of capture", ylim = c(0, 1))
## for(i in 1:k){
## 	lines(c(years[i], years[i]), c(an.capt.q[2,i], an.capt.q[3,i]) )
## }

## #Figura español:
## plot(years, an.capt.q[1,], pch = 15, xlab = "Año(verano)", ylab = "Probabilidad  de captura", ylim = c(0, 1.0))
## for(i in 1:k){
## 	lines(c(years[i], years[i]), c(an.capt.q[2,i], an.capt.q[3,i]) )
## }

## write.table(an.capt.q, "pr_captura.csv", sep = ',')

## #POST-HOC RELATIONSHIP OF EFFORT AND PR. CAPTURE
## plot(raw.effort, an.capt.q[1,], pch = 15, xlab = "Effort (days)", ylab = "Pr. of capture", ylim = c(0, 1))
## for(i in 1:k){
## 	lines(c(raw.effort[i], raw.effort[i]), c(an.capt.q[2,i], an.capt.q[3,i]) )
## }


## #PR OF RECAPTURE AS FCN OF EFFORT - when using effort models
## new.effort<-seq(-2,2, by = 0.1)
## #chain values for gamma0 and gamma1
## gam0.ind<-which(varnames(mod1) == "gamma[1]")
## the.gam0<-append( mod1[[3]][,gam0.ind], append(mod1[[1]][,gam0.ind], mod1[[2]][,gam0.ind]) )

## gam1.ind<-which(varnames(mod1) == "gamma[2]")
## the.gam1<-append( mod1[[3]][,gam1.ind], append(mod1[[1]][,gam1.ind], mod1[[2]][,gam1.ind]) )

## new.effort.mat<-matrix(nrow=length(the.gam0), ncol = length(new.effort))
## new.effort.mat.q<-matrix(nrow = 3, ncol = length(new.effort))
## for(j in 1:length(new.effort)){
## 	new.effort.mat[,j]<- inv.logit(the.gam0 + the.gam1*new.effort[j])
## 	new.effort.mat.q[,j]<-quantile(new.effort.mat[,j], probs = c(0.5, 0.025, 0.975) )
## }

## plot(new.effort*effort.sd + effort.bar, new.effort.mat.q[1,], type = 'l', ylim = c(0,1), xlab = "Effort (days)", ylab = "Probability of capture")
## lines(new.effort*effort.sd + effort.bar, new.effort.mat.q[2,], lty = 2)
## lines(new.effort*effort.sd + effort.bar, new.effort.mat.q[3,], lty = 2)
## points(y=rep(0.01,length(raw.effort)) , x=raw.effort, pch = 17)

## #look at relationship between etap and effort
## some.etap<-summary(mod1[,which(varnames(mod1) =="etap[1]"):which(varnames(mod1) == "etap[9]"), ], quantiles = c(0.5, 0.025, 0.975) )[[2]]

## plot(effort, some.etap[,1], pch = 15)












