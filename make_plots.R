### Make some basic plots of the different model versions  to comare them

## The fits arranged into a list for easier processing below
posts <- list(post.log1, post.vb1, post.vb_surv, post.vb_surv_tvzeta)
model.names <- c('log1', 'vb1', 'vb_surv', 'vb_surv_tvzeta')

## Look at time series of abundance
out <- ldply(1:4, function(i)
  cbind(model= model.names[i], years=years,
        ddply(melt(posts[[i]][,paste("N", 1:nyrs, sep='_')], id.vars=NULL), .(variable), summarize,
              lwr=quantile(value, .025),
              med=median(value),
              upr=quantile(value, .975))))
g <- ggplot(out, aes(x=years)) + geom_ribbon(aes(ymin=lwr, ymax=upr), fill=gray(.6)) +
  geom_line(aes(y=med), lwd=2) + facet_wrap('model')
ggsave('plots/abundances.png', g, width=ggwidth, height=ggheight)

## Abundance in last year
out <- ldply(1:4, function(i)
  data.frame(model= model.names[i], terminalN=posts[[i]][,paste0('N_',nyrs)]))
g <- ggplot(out, aes(x=model,y=terminalN)) + geom_violin()
ggsave('plots/terminal_abundances.png', g, width=ggwidth, height=ggheight)

## Survival probability
out <- ldply(1:4, function(i)
  data.frame(model= model.names[i], survival=posts[[i]][,'surv']))
g <- ggplot(out, aes(x=model,y=survival)) + geom_violin()
ggsave('plots/survival.png', g, width=ggwidth, height=ggheight)
