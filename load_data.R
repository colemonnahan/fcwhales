### Read in raw data from files
message("Loading real sightings data...")


if(!file.exists("data/capture_histories_2017_11_28.csv"))
  stop("Data file does not exist")
captures <- read.csv("data/capture_histories_2017_11_28.csv") # capture histories
n.obs <- dim(captures)[1] # number of unique whales
## Make sure each year has 8 columns and each month has 15 years
tt <- names(captures)
year.check <- any(nseas !=as.numeric(table(as.numeric(substr(tt,2,5)))))
if(year.check) stop("Wrong number of seasons for years in capture data")
month.check <- any(nyrs != as.numeric(table(substr(tt,7,9)) ))
if(month.check) stop("Wrong number of years in a season in capture data")
if(n.timesteps != dim(captures)[2])
  stop("wrong number of columns in capture history data file")

### Process data for modeling and plotting
## Add extra rows to make space for the latent whales
M <- n.obs+50
post <- matrix(nrow=M, ncol=n.timesteps, data=0)
post[1:n.obs, 1:n.timesteps] <- as.matrix(captures)
## Convert capture histories into array form
x <- array(data=NA, dim=c(M, nyrs, nseas))
dimnames(x) <- list("ID"=1:dim(x)[1], "year"=years, "month"=months)
for(i in 0:(nyrs-1) ){
	x[ ,i+1 , ] <- post[, (i*nseas + 1):(i*8 + nseas)]
}
if(sum(is.na(x)) > 0)
  stop("Issue converting captures to array -- NAs present")

message("Loading effort data...")
## Effort:year x month. Each cells is the number of days of effort.
pre.effort <- read.table("data/effort_2017_11_28.csv", header=F, sep=',')
## transpose so month x year
effort <- as.matrix(pre.effort)
dimnames(effort) <- list("year"=years, "month"=months)
effort.long <- melt(effort, value.name='effort')


message("Checking for data issues...")
captures.long <- melt(x)
captures.long <- ddply(captures.long, .(month, year), summarize, unique.whales=sum(value))
## Use this to check back with original data if needed
## captures.long <- captures.long[order(captures.long$year, captures.long$month),]
obs.effort <- merge(effort.long, captures.long, by=c("year", "month"))
effort.check <- droplevels(subset(obs.effort, effort == 0 & unique.whales > 0))
if(nrow(effort.check)>0){
  print(effort.check)
  warning("These months have 0 effort but observed whales")
}

## This is whether each whale was seen in a year
seen.years <- apply(x, c(1,2), max)
## Now get the first year seen
temp <- as.numeric(apply(seen.years, 1, function(ii) which(ii==1)[1]))
## Drop the NAs and convert to year.
first.year.seen <- years[temp[1:nrow(captures)]]
## Trick to get table to count 0 whales when year is missing
first.year.seen <- factor(first.year.seen, levels=years)
x1 <- as.numeric(table(first.year.seen))
x2 <- as.numeric(apply(apply(x, c(1,2) , max), 2, sum)) # unique whales by year
e <- as.numeric(apply(effort, 1, sum))
whales.table <- data.frame(year=years+1, effort=e, unique.whales=x2, new.whales=x1)

### No longer using this.
## ## start months and end months for each visit
## non.visit <- read.table("data/not_visited15.csv", header=T, row.names=1, sep=',')
## if(nseas != dim(non.visit)[1])
##   stop("wrong number of rows in non_visit data file")
## NV <- as.matrix(non.visit)
## dimnames(NV) <- list("month"=months, "year"=years)

## Prep the latent whale JAGS inputs
w <- rep(NA, length=M)
w[1:n.obs] <- 1
## Real data for all models
realdata <- list('w'=w, 'M'=M, 'nyrs'=nyrs, 'nseas'=nseas, 'x'= x,
                 'effort'=effort)

message("Dropping years 2002, 2003, and 2007 from data")
## We decided to drop the first two years since the effort is really
## different from other years. Earlier tests showed this had little impact
## on key outputs like survival, but should make the availability estimates
## better.
realdata$nyrs <- realdata$nyrs-2
realdata$x <- realdata$x[,-(1:2),]
realdata$effort <- realdata$effort[-(1:2),]
## We also decided to take out the year 2007/2008 since the photos were
## made by someone else and dont really fit with the realdataa collection
## scheme.
realdata$effort["2007",] <- 0
realdata$x[,"2007",] <- 0
nyrs <- nyrs-2
years <- 2003+1:nyrs
