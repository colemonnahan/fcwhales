

## This is a copy of the real effort data which are not confidential.
effort.sim <-
  structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 2L,
              0L, 0L, 0L, 3L, 9L, 3L, 3L, 0L, 2L, 4L, 0L, 6L, 0L, 4L, 7L, 9L,
              12L, 0L, 4L, 4L, 0L, 0L, 0L, 5L, 11L, 6L, 6L, 0L, 6L, 7L, 8L,
              14L, 0L, 6L, 13L, 7L, 6L, 0L, 7L, 2L, 7L, 8L, 0L, 5L, 8L, 9L,
              14L, 3L, 3L, 7L, 0L, 4L, 3L, 4L, 1L, 3L, 8L, 0L, 0L, 5L, 11L,
              14L, 0L, 0L, 2L, 5L, 4L, 3L, 3L, 0L, 4L, 4L, 0L, 0L, 7L, 8L,
              6L, 0L, 0L, 0L, 5L, 3L, 0L, 1L, 1L, 3L, 2L, 3L, 0L, 0L, 0L, 4L,
              0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
              0L, 0L, 0L), .Dim = c(15L, 8L),
            .Dimnames = structure(list(
              year = c("2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                       "2013", "2014", "2015", "2016"),
              month = c("Nov", "Dec", "Jan", "Feb", "Mar",
                        "Apr", "May", "Jun")),
              .Names = c("year", "month")))

pa1 <- rnorm(nyrs, 2,.3)
pa1[5] <- .5 # one really low year to see if it can recover this
pars <- list(piab=0.2, surv=.9, pa1=inv.logit(pa1),
             pa2=c(.05,.1, .2, .4, .5, .3, .2, .08),
             k=.2)
simdata <- simulate_data(N=125, pars=pars, effort=effort.sim,
                         seed=1523, nyrs=nyrs, nseas=nseas)
dimnames(simdata$dat$x) <-
  list(ID=as.character(1:nrow(simdata$dat$x)), year=dimnames(effort.sim)[[1]], month=dimnames(effort.sim)[[2]])
plot.simdat(simdata)
