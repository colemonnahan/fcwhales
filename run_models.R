## This file runs models on the real data and saves the fits to file for
## later processing.

## We decided to drop the first two year sicne the effort is really
## different from other years. Earlier tests showed this had little impact
## on key outputs like survival, but should make the availability estimates
## better.
message("Dropping years 2002, 2003, and 2007 from data")
nrem <- 2
dat2 <- dat
dat2$nyrs <- dat2$nyrs-nrem
dat2$x <- dat2$x[,-(1:nrem),]
dat2$effort <- dat2$effort[-(1:nrem),]
## We also decided to take out the year 2007/2008 since the photos were
## made by someone else and dont really fit with the data collection
## scheme.
dat2$effort["2007",] <- 0
dat2$x[,"2007",] <- 0

## MCMC settings; th=thing; ni= iterations (total); nc= chains (parallel)
th <- 1; ni <- 100*th; nc <- 3; (ni/2/th*nc) # quick laptop runs
th <- 100; ni <- 500*th; nc <- 4; (ni/2/th*nc) # full runs

t1 <- Sys.time()
fit.vb1 <- run_vb(dat2, ni, th, nc)
saveRDS(fit.vb1, file='results/fit.vb1.RDS')
message("Model vb1:"); print(Sys.time()-t1)

t1 <- Sys.time()
fit.vb_surv <- run_vb(dat2, ni, th, nc, surv='informative')
saveRDS(fit.vb_surv, file='results/fit.vb_surv.RDS')
message("Model vb_surv:"); print(Sys.time()-t1)

t1 <- Sys.time()
fit.log1 <- run_logistic(dat2, ni, th, nc)
saveRDS(fit.log1, file='results/fit.log1.RDS')
message("Model log1:"); print(Sys.time()-t1)

## ## Same models with full data. Despite the poor naming convention. Sorry!
## t1 <- Sys.time()
## fit.vb2 <- run_vb(dat, ni, th, nc)
## saveRDS(fit.vb2, file='results/fit.vb2.RDS')
## message("Model log2:"); print(Sys.time()-t1)
## t1 <- Sys.time()
## fit.log2 <- run_logistic(dat, ni, th, nc)
## saveRDS(fit.log2, file='results/fit.log2.RDS')
## message("Model log2:"); print(Sys.time()-t1)



