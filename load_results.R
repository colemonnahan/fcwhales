## This file loads in the saved results and processes them where needed

fit.log1 <- readRDS(file='results/fit.log1.RDS')
fit.vb1 <- readRDS(file='results/fit.vb1.RDS')
fit.vb_surv <- readRDS(file='results/fit.vb_surv.RDS')
fit.vb_surv_tvzeta <- readRDS(file='results/fit.vb_surv_tvzeta.RDS')

## Some convenient data.frames
post.vb1 <- convert_to_df(fit.vb1, 'vb1')
post.log1 <- convert_to_df(fit.log1, 'log1')
post.vb_surv <- convert_to_df(fit.vb_surv, 'vb1surv')
post.vb_surv_tvzeta <- convert_to_df(fit.vb_surv_tvzeta, 'vb_surv_tvzeta')
