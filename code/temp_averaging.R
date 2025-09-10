library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')
library('paletteer')

k <- 6.791884
setwd("/Users/cruzloya/git/zikvTPCnla/")

data.constT <- read.csv("./data/processed_data/constantT.csv")
data.fluctmax <- read.csv("./data/processed_data/fluctmax.csv")
data.fluctmin <- read.csv("./data/processed_data/fluctmin.csv")

temp.data <- read.csv("./data/temp_program/temp_program.csv")

mean.temps <- unique(temp.data$mean.temp)
cols <- rev(paletteer_c("ggthemes::Orange-Blue Diverging",
                    length(mean.temps))) 

get.temps <- function(times, known.times, known.temps) {
  f <- approxfun(c(known.times, 24), c(known.temps, known.temps[1]))
  t <- times - (24 * times %/% 24)
  return(f(t))
}

## Plot temperature history.
par(mfrow=c(1,1))
## Constant
times <- seq(0, 72, 0.1)
plot("", xlim=c(0, 72), ylim=c(10, 40), xlab="time [hour]",
     ylab="temperature [°C]")
for(i in 1:length(mean.temps)) {
  temp <- mean.temps[i]
  lines(times, rep(temp, length(times)), col=cols[i], lwd=2)
}


## 5 DTR
times <- seq(0, 72, 0.1)
plot("", xlim=c(0, 72), ylim=c(10, 40), xlab="time [hour]",
     ylab="temperature [°C]")
for(i in 1:length(mean.temps)) {
  temp <- mean.temps[i]
  df <- subset(temp.data, (temp.data$mean.temp == temp) & 
                          (temp.data$DTR == 5))
  temps <- get.temps(times, df$time, df$temperature)
  lines(times, temps, col=cols[i], lwd=2)
}

## 10 DTR
times <- seq(0, 72, 0.1)
plot("", xlim=c(0, 72), ylim=c(10, 40), xlab="time [hour]",
     ylab="temperature [°C]")
for(i in 1:length(mean.temps)) {
  temp <- mean.temps[i]
  df <- subset(temp.data, (temp.data$mean.temp == temp) & 
                 (temp.data$DTR == 10))
  temps <- get.temps(times, df$time, df$temperature)
  lines(times, temps, col=cols[i], lwd=2)
}

constantT <- readRDS("./MCMC_chains/constant_temp.RDS")
fluctmin <- readRDS("./MCMC_chains/fluctmin_temp.RDS")
fluctmax <- readRDS("./MCMC_chains/fluctmax_temp.RDS")

# FlexTPC model for thermal performance curves.
flexTPC <- function(temp, Tmin, Tmax, rmax, alpha, beta) {
  s <- alpha * (1 - alpha) / beta^2
  result <- rep(0, length(temp))
  Tidx = (temp > Tmin) & (temp < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (temp[Tidx] - Tmin) / alpha) 
                                  + (1 - alpha) * log( (Tmax - temp[Tidx]) / (1 - alpha)) 
                                  - log(Tmax - Tmin)) ) 
  return(result)
}


temps <-seq(0, 40, 0.1)
plot(temps, flexTPC(temps, 10, 35, 1, 0.8, 0.2), type='l')

## Transformation and inverse transformation functions.
f <- function(y, k=6.791884) {
  return(log10(1 + y / 10^k))
}

f.inv <- function(z, k=6.791884) {
  return(10^(z+k) - 10^k)
}


constantT.chains <- MCMCchains(constantT, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
fluctmin.chains <- MCMCchains(fluctmin, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))
fluctmax.chains <- MCMCchains(fluctmax, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))


temps <- seq(10, 35, 0.1)
constantT.curves <- apply(constantT.chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
fluctmin.curves <- apply(fluctmin.chains, 1, 
                         function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))
fluctmax.curves <- apply(fluctmax.chains, 1, 
                         function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))


## 5 DTR / fluctmin
times <- seq(0, 24, 0.1)
# 18C data
df <- subset(temp.data, (temp.data$mean.temp == 18) & 
               (temp.data$DTR == 5))
temps.18 <- get.temps(times, df$time, df$temperature)

mean.temps
nla.fluctmin.preds <- matrix(nrow=200000, ncol=length(mean.temps))
for(i in 1:length(mean.temps)) {
  temp <- mean.temps[i]
  temps <- temps.18 + (temp - 18)
  pred.posterior.mu <- apply(constantT.chains, 1, function(x) 
    f.inv(flexTPC(temps, x[1], x[2], x[3], x[4], x[5])))
  nla.fluctmin.preds[, i] <- f(apply(pred.posterior.mu, 2, mean))
}


plot(data.fluctmin$temp, log10(1 + data.fluctmin$vir.rate / 10^k), pch=20,
     ylim=c(0, 4.5), xlab='temperature [°C]',
     ylab='log(1 + 10^-k * growth rate)',
     main='2.5°C fluctuation')

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(constantT.curves, 1, mean)
CI  <- apply(constantT.curves, 1, quantile, c(0.025, 0.975))

temps <- seq(10, 35, 0.1)

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="steelblue", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("steelblue", 0.2), lty=0)

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(fluctmin.curves, 1, mean)
CI  <- apply(fluctmin.curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="darkgreen", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("darkgreen", 0.2), lty=0)

lines(mean.temps, apply(nla.fluctmin.preds, 2, mean), col='purple', lwd=1.5)
polygon(c(mean.temps, rev(mean.temps)), c(apply(nla.fluctmin.preds, 2, quantile, 0.025), 
                                          rev(apply(nla.fluctmin.preds, 2, quantile, 0.975))), 
        col=alpha("purple", 0.2), lty=0)


## 10 DTR / fluctmax
times <- seq(0, 24, 0.1)
# 18C data
df <- subset(temp.data, (temp.data$mean.temp == 18) & 
               (temp.data$DTR == 10))
temps.18 <- get.temps(times, df$time, df$temperature)

mean.temps
nla.fluctmax.preds <- matrix(nrow=200000, ncol=length(mean.temps))
for(i in 1:length(mean.temps)) {
  temp <- mean.temps[i]
  temps <- temps.18 + (temp - 18)
  pred.posterior.mu <- apply(constantT.chains, 1, function(x) 
    f.inv(flexTPC(temps, x[1], x[2], x[3], x[4], x[5])))
  nla.fluctmax.preds[, i] <- f(apply(pred.posterior.mu, 2, mean))
}


plot(data.fluctmax$temp, log10(1 + data.fluctmax$vir.rate / 10^k), pch=20,
     ylim=c(0, 4.5), xlab='temperature [°C]',
     ylab='log(1 + 10^-6 * growth rate)', 
     main='5°C fluctuation')

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(constantT.curves, 1, mean)
CI  <- apply(constantT.curves, 1, quantile, c(0.025, 0.975))

temps <- seq(10, 35, 0.1)

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="steelblue", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("steelblue", 0.2), lty=0)

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(fluctmax.curves, 1, mean)
CI  <- apply(fluctmax.curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="darkgreen", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("darkgreen", 0.2), lty=0)

lines(mean.temps, apply(nla.fluctmax.preds, 2, mean), col='purple', lwd=1.5)
polygon(c(mean.temps, rev(mean.temps)), c(apply(nla.fluctmax.preds, 2, quantile, 0.025), 
                                          rev(apply(nla.fluctmax.preds, 2, quantile, 0.975))), 
        col=alpha("purple", 0.2), lty=0)



# Fluctmin predictions. 
constantT.preds <- apply(constantT.chains, 1, 
                          function(x) flexTPC(data.fluctmin$temp, x[1], x[2], x[3], x[4], x[5]))
fluctmin.preds <- apply(fluctmin.chains, 1, 
                         function(x) flexTPC(data.fluctmin$temp, x[1], x[2], x[3], x[4], x[5]))

y_obs <- log10(1 + data.fluctmin$vir.rate / 10^k)
str(y_obs)

constantT.preds
data.fluctmin
constantT.err <- apply((constantT.preds - y_obs)^2, 2, sum)
fluctmin.err <- apply((fluctmin.preds - y_obs)^2, 2, sum)
nla.fluctmin.err <- apply((t(cbind(nla.fluctmin.preds[, c(1,2,3,6,7,8)], 
                                   nla.fluctmin.preds[, c(1,2,3,6,7,8)], 
                                   nla.fluctmin.preds[, c(1,2,3,6,7,8)]))
                           - y_obs)^2, 2, sum)

par(mfrow=c(1, 2))
plot("", xlim=c(1, 3), ylim=c(0, 15), ylab="mean squared error", xlab="")

rel.err.constantT <- constantT.err / fluctmin.err
rel.err.nla <- nla.fluctmin.err / fluctmin.err
points(c(1, 2, 3), c(mean(constantT.err), mean(nla.fluctmin.err),
                     mean(fluctmin.err)), pch=20)
points(c(1, 2, 3), c(quantile(constantT.err, 0.025), 
                     quantile(nla.fluctmin.err, 0.025),
                     quantile(fluctmin.err, 0.025)), pch=20)
points(c(1, 2, 3), c(quantile(constantT.err, 0.975), 
                     quantile(nla.fluctmin.err, 0.975),
                     quantile(fluctmin.err, 0.975)), pch=20)
#abline(h=1, lty=2)


str(constantT.preds)  
  
mean((constantT.preds - y_obs)^2)
mean((t(cbind(nla.fluctmin.preds[, c(1,2,3,6,7,8)],
              nla.fluctmin.preds[, c(1,2,3,6,7,8)], 
              nla.fluctmin.preds[, c(1,2,3,6,7,8)]))
      - y_obs)^2)
mean((fluctmin.preds - y_obs)^2)


# Fluctmax predictions. 
constantT.preds <- apply(constantT.chains, 1, 
                         function(x) flexTPC(data.fluctmax$temp, x[1], x[2], x[3], x[4], x[5]))
fluctmax.preds <- apply(fluctmax.chains, 1, 
                        function(x) flexTPC(data.fluctmax$temp, x[1], x[2], x[3], x[4], x[5]))

y_obs <- log10(1 + data.fluctmax$vir.rate / 10^k)
str(y_obs)

mean((constantT.preds - y_obs)^2)
mean((t(cbind(nla.fluctmax.preds[, c(1,2,3,6,7,8)], 
              nla.fluctmax.preds[, c(1,2,3,6,7,8)], 
              nla.fluctmax.preds[, c(1,2,3,6,7,8)]))
      - y_obs)^2)
mean((fluctmax.preds - y_obs)^2)






####
par(mfrow=c(1,3))

plot(data.constT$temp, log10(1 + data.constT$vir.rate), pch=20,
     ylim=c(0, 4.5), xlab='temperature [°C]', xlim=c(12, 33),
     ylab='log(1 + 10^-6 * growth rate)',
     main='constant temperature')

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(constantT.curves, 1, mean)
CI  <- apply(constantT.curves, 1, quantile, c(0.025, 0.975))

temps <- seq(10, 35, 0.1)

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="steelblue", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("steelblue", 0.2), lty=0)


plot(data.fluctmin$temp, log10(1 + data.fluctmin$vir.rate), pch=20,
     ylim=c(0, 4.5), xlab='temperature [°C]', xlim=c(12, 33),
     ylab='log(1 + 10^-6 * growth rate)',
     main='2.5°C fluctuation')

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(fluctmin.curves, 1, mean)
CI  <- apply(fluctmin.curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="darkgreen", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("darkgreen", 0.2), lty=0)


plot(data.fluctmax$temp, log10(1 + data.fluctmax$vir.rate), pch=20,
     ylim=c(0, 4.5), xlab='temperature [°C]', xlim=c(12, 33),
     ylab='log(1 + 10^-6 * growth rate)',
     main='5°C fluctuation')



# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(fluctmax.curves, 1, mean)
CI  <- apply(fluctmax.curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="darkgreen", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("darkgreen", 0.2), lty=0)


