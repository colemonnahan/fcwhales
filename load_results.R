## This file loads in the saved results and processes them where needed

fit.log1 <- readRDS(file='results/fit.log1.RDS')
fit.log2 <- readRDS(file='results/fit.log2.RDS')
fit.vb1 <- readRDS(file='results/fit.vb1.RDS')
fit.vb_surv <- readRDS(file='results/fit.vb_surv.RDS')
fit.vb2 <- readRDS(file='results/fit.vb2.RDS')

## Some convenient data.frames
post.vb1 <- convert_to_df(fit.vb1, 'vb1')
post.vb2 <- convert_to_df(fit.vb2, 'vb2')
post.log1 <- convert_to_df(fit.log1, 'log1')
post.vb_surv <- convert_to_df(fit.vb_surv, 'vb1surv')
post.log2 <- convert_to_df(fit.log2, 'log2')
