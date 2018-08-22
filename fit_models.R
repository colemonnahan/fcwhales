## This file runs models on the real data and saves the fits to file for
## later processing.

## MCMC settings; th=thing; ni= iterations (total); nc= chains (parallel)
th <- 1; ni <- 100*th; nc <- 3; (ni/2/th*nc) # quick laptop runs
## th <- 100; ni <- 500*th; nc <- 4; (ni/2/th*nc) # full runs

t1 <- Sys.time()
fit.vb <- run_vb(dat, ni, th, nc)
saveRDS(fit.vb, file='results/fit.vb.RDS')
message("Model vb:"); print(Sys.time()-t1)

t1 <- Sys.time()
fit.vb_surv <- run_vb(dat, ni, th, nc, surv='informative')
saveRDS(fit.vb_surv, file='results/fit.vb_surv.RDS')
message("Model vb_surv:"); print(Sys.time()-t1)

t1 <- Sys.time()
fit.vb_surv_tvzeta <-
  run_vb(dat, ni, th, nc, surv='informative', tvzeta=TRUE)
saveRDS(fit.vb_surv_tvzeta, file='results/fit.vb_surv_tvzeta.RDS')
message("Model vb_surv_tvzeta:"); print(Sys.time()-t1)

t1 <- Sys.time()
fit.log <- run_logistic(dat, ni, th, nc)
saveRDS(fit.log, file='results/fit.log.RDS')
message("Model log:"); print(Sys.time()-t1)




