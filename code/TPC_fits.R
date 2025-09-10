library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

setwd("/Users/cruzloya/git/zikvTPCnla/")

# MLE estimate for transformation
k <- 6.791884

data.constT <- read.csv("./data/processed_data/constantT.csv")
data.fluctmax <- read.csv("./data/processed_data/fluctmax.csv")
data.fluctmin <- read.csv("./data/processed_data/fluctmin.csv")


set.seed(42)

data.constT$vir.rate
data.fluctmax$vir.rate
par(mfrow=c(1,3))
plot(data.constT$temp, data.constT$vir.rate * 10^-9, xlab="temperature",
     ylab="viruses / hour * 10^-9", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 15), pch=20, col="steelblue")
plot(data.fluctmin$temp, data.fluctmin$vir.rate* 10^-9, xlab="temperature", 
     ylab="viruses / hour * 10^-9", 
     main="2.5C fluctuation", xlim=c(15, 38), ylim=c(0, 15), pch=20, col="steelblue")
plot(data.fluctmax$temp, data.fluctmax$vir.rate* 10^-9, xlab="temperature", 
     ylab="viruses / hour * 10^-9", 
     main="5C fluctuation", xlim=c(15, 38), ylim=c(0, 15), pch=20, col="steelblue")

par(mfrow=c(1,3))
plot(data.constT$temp, log10(1 + data.constT$vir.rate / 10^k), xlab="temperature",
     ylab="log(1 + viral titer / hour * 10^-6)", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")
plot(data.fluctmin$temp, log10(1 + data.fluctmin$vir.rate / 10^k), xlab="temperature", 
     ylab="log(1 + viral titer / hour * 10^-6)", 
     main="2.5C fluctuation", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")
plot(data.fluctmax$temp, log10(1 + data.fluctmax$vir.rate / 10^k), xlab="temperature", 
     ylab="log(1 + viral titer / hour * 10^-6)", 
     main="5C fluctuation", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")


# FlexTPC model for thermal performance curves.
flexTPC <- function(T, Tmin, Tmax, rmax, alpha, beta) {
  s <- alpha * (1 - alpha) / beta^2
  result <- rep(0, length(T))
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + (1 - alpha) * log( (Tmax - T[Tidx]) / (1 - alpha))
                                  - log(Tmax - Tmin)) ) 
  return(result)
}

## Defines JAGS model.
sink("flexTPC_norm.txt")
cat("
    model{
    
    ## Priors (please make sure to change these priors so they make sense for your system!)
    Tmin ~ dnorm(15, 1/5^2) # Minimum temperature (Celsius). Approx. prior 95% CI: [5, 25] 
    Tmax ~ dnorm(35, 1/5^2) # Maximum temperature (Celsius). Approx. prior 95% CI: [25, 45]
    rmax ~ dunif(0, 10)  
    alpha ~ dunif(0, 1) # Relative position of Topt relative to Tmin and Tmax.
                        # Uniform prior gives equal prior weight to curves of any skewness.
                        # Can also use e.g. beta prior to give more prior probability to left-skewed, symmetric
                        # or right-skewed curves depending on parameters of prior if it makes sense for the trait.
    beta ~ dgamma(0.3^2/0.3^2, 0.3/0.3^2) # Gamma prior for upper thermal breath. Values
                                          # around 0.2-0.4 correspond to common TPC shapes
                                          # as are described by e.g. the Briere1 or quadratic
                                          # models. 95% prior CI: [0.008, 1.107]
    sigma ~ dunif(0, 2) # Standard deviation for the data points around the fitted TPC. 
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    B <- beta * (Tmax - Tmin)

    ## Likelihood
    for(i in 1:N.obs){
      # flexTPC model.
      mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                            + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                            - log(Tmax - Tmin)) ) 
      # Normal likelihood with constTant variance.
      
      y[i] ~ dnorm(mu[i], 1 / sigma^2)T(0, )
    }
    
    } # close model
    ",fill=TRUE)
sink()

## Defines JAGS model.
sink("flexTPC_norm_var_SD.txt")
cat("
    model{
    
    ## Priors (please make sure to change these priors so they make sense for your system!)
    Tmin ~ dnorm(15, 1/5^2) # Minimum temperature (Celsius). Approx. prior 95% CI: [5, 25] 
    Tmax ~ dnorm(35, 1/5^2) # Maximum temperature (Celsius). Approx. prior 95% CI: [25, 45]
    rmax ~ dunif(0, 10)  
    alpha ~ dunif(0, 1) # Relative position of Topt relative to Tmin and Tmax.
                        # Uniform prior gives equal prior weight to curves of any skewness.
                        # Can also use e.g. beta prior to give more prior probability to left-skewed, symmetric
                        # or right-skewed curves depending on parameters of prior if it makes sense for the trait.
    beta ~ dgamma(0.3^2/0.3^2, 0.3/0.3^2) # Gamma prior for upper thermal breath. Values
                                          # around 0.2-0.4 correspond to common TPC shapes
                                          # as are described by e.g. the Briere1 or quadratic
                                          # models. 95% prior CI: [0.008, 1.107]
    sigma0 ~ dunif(0, 1) # Standard deviation for the data points around the fitted TPC. 
                        # Upper limit of 1 for this problem, but please change this
                        # to make sense for the data you're fitting.
    sigma_beta ~ dnorm(0, 1)
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    B <- beta * (Tmax - Tmin)

    ## Likelihood
    for(i in 1:N.obs){
      # flexTPC model.
      mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                            + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                            - log(Tmax - Tmin)) ) 
      # Normal likelihood with constTant variance.
      
      sigma[i] <- sigma0 * exp(sigma_beta * mu[i])
      y[i] ~ dnorm(mu[i], 1 / sigma[i]^2)T(0, )
    }
    
    } # close model
    ",fill=TRUE)
sink()

sink("flexTPC_gamma.txt")
cat("
    model{
    
    ## Priors (please make sure to change these priors so they make sense for your system!)
    Tmin ~ dnorm(15, 1/5^2) # Minimum temperature (Celsius). Approx. prior 95% CI: [5, 25] 
    Tmax ~ dnorm(35, 1/5^2) # Maximum temperature (Celsius). Approx. prior 95% CI: [25, 45]
    rmax ~ dexp(1 / 10)  
    alpha ~ dunif(0, 1) # Relative position of Topt relative to Tmin and Tmax.
                        # Uniform prior gives equal prior weight to curves of any skewness.
                        # Can also use e.g. beta prior to give more prior probability to left-skewed, symmetric
                        # or right-skewed curves depending on parameters of prior if it makes sense for the trait.
    beta ~ dgamma(0.3^2/0.3^2, 0.3/0.3^2) # Gamma prior for upper thermal breath. Values
                                          # around 0.2-0.4 correspond to common TPC shapes
                                          # as are described by e.g. the Briere1 or quadratic
                                          # models. 95% prior CI: [0.008, 1.107]
    sigma0 ~ dunif(0, 10) # Standard deviation for the data points around the fitted TPC. 
                        # Upper limit of 1 for this problem, but please change this
                        # to make sense for the data you're fitting.
    sigma_beta ~ dnorm(0, 1)
    
    # Derived quantities
    s <- alpha * (1 - alpha) / beta^2
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    B <- beta * (Tmax - Tmin)

    ## Likelihood
    for(i in 1:N.obs){
      # flexTPC model.
      mu[i] <- (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (alpha * log( max(temp[i] - Tmin, 10^-20) / alpha) 
                            + (1 - alpha) * log( max(Tmax - temp[i], 10^-20) / (1 - alpha))
                            - log(Tmax - Tmin)) ) 
      # Normal likelihood with constTant variance.
      
      sigma[i] <- sigma0 * exp(sigma_beta * mu[i])
      y[i] ~ dgamma(mu[i]^2 / sigma[i]^2, mu[i] / sigma[i]^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### Set MCMC Settings
ni <- 500000 # number of iterations in each chain
nb <- 100000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 4 # number of chains

##### constTant data

##### Organize Data for JAGS
y <- log10(1 + data.constT$vir.rate / 10^k)
temp <- data.constT$temp
N.obs <- length(y)

##### Bundle Data
jag.data<-list(y=y, temp = temp, N.obs=N.obs)

### Initial values for MCMC chains (chosen randomly from specified 
### distributions). Please remember to change these to something suitable
### for the data you're fitting.

inits<-function(){list(
  Tmin = runif(1, min=5, max=20),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=3.0, max=8.0),
  alpha = runif(1, min=0, max=1),
  beta = runif(1, 0.1, 0.5),
  sigma = runif(1, min=0, max=1))}


##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "B",
                "sigma")

jags.out.constT <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                 model.file="flexTPC_norm.txt", n.thin=nt, n.chains=nc, 
                 n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
jags.out.constT

mcmcplot(jags.out.constT)


# Plot data.
par(mfrow=c(1,3))
plot(data.constT$temp, log10(1 + data.constT$vir.rate / 10^k), xlab="temperature",
     ylab="log(1 + viral titer / hour * 10^k)", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out.constT, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="orange", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("orange", 0.2), lty=0)


### Fluctmin

##### Organize Data for JAGS
y <- log10(1 + data.fluctmin$vir.rate / 10^k)
temp <- data.fluctmin$temp
N.obs <- length(y)

##### Bundle Data
jag.data<-list(y=y, temp = temp, N.obs=N.obs)

### Initial values for MCMC chains (chosen randomly from specified 
### distributions). Please remember to change these to something suitable
### for the data you're fitting.

inits<-function(){list(
  Tmin = runif(1, min=5, max=20),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=3.0, max=8.0),
  alpha = runif(1, min=0, max=1),
  beta = runif(1, 0.1, 0.5),
  sigma = runif(1, min=0, max=1))}


##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "B",
                "sigma")

jags.out.fluctmin <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                 model.file="flexTPC_norm.txt", n.thin=nt, n.chains=nc, 
                 n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
jags.out.fluctmin

mcmcplot(jags.out.fluctmin)

# Plot data.
plot(data.fluctmin$temp, log10(1 + data.fluctmin$vir.rate / 10^k), xlab="temperature", 
     ylab="log(1 + viral titer / hour * 10^-k)", 
     main="2.5C fluctuation", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out.fluctmin, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="orange", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("orange", 0.2), lty=0)

### Fluctmax

##### Organize Data for JAGS
y <- log10(1 + data.fluctmax$vir.rate / 10^k)
temp <- data.fluctmax$temp
N.obs <- length(y)

##### Bundle Data
jag.data<-list(y=y, temp = temp, N.obs=N.obs)

### Initial values for MCMC chains (chosen randomly from specified 
### distributions). Please remember to change these to something suitable
### for the data you're fitting.

inits<-function(){list(
  Tmin = runif(1, min=5, max=20),
  Tmax = runif(1, min=30, max=40),
  rmax = runif(1, min=3.0, max=8.0),
  alpha = runif(1, min=0, max=1),
  beta = runif(1, 0.1, 0.5),
  sigma = runif(1, min=0, max=1))}


##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt", "B",
                "sigma")

jags.out.fluctmax <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                          model.file="flexTPC_norm.txt", n.thin=nt, n.chains=nc, 
                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
jags.out.fluctmax

mcmcplot(jags.out.fluctmax)

# Plot data.
plot(data.fluctmax$temp, log10(1 + data.fluctmax$vir.rate / 10^k), xlab="temperature", 
     ylab="log(1 + viral titer / hour * 10^-6)", 
     main="5C fluctuation", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out.fluctmax, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="orange", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("orange", 0.2), lty=0)

View(jags.out.constT$BUGSoutput$summary)
View(jags.out.fluctmin$BUGSoutput$summary)
View(jags.out.fluctmax$BUGSoutput$summary)




# Multiple TPCs in same plot.
par(mfrow=c(1,1))
plot("", xlab="temperature", 
     ylab="log(1 + viral titer / hour * 10^-k)", 
     main="Constant vs fluctuating", xlim=c(15, 40), ylim=c(0, 4))

### Constant

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out.constT, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="steelblue", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("steelblue", 0.2), lty=0)

### 2.5C

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out.fluctmin, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="orange", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("orange", 0.2), lty=0)

# 5C

# Extract MCMC chains for flexTPC parameters.
chains <- MCMCchains(jags.out.fluctmax, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))

# Calculate TPCs for each posterior sample.
temps <- seq(0, 40, 0.1)
curves <- apply(chains, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Calculate posterior mean of TPCs at every temperature and 95% credible interval.
meancurve <- apply(curves, 1, mean)
CI  <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot posterior mean and 95% credible interval.
lines(temps, meancurve, col="darkgreen", lwd=1.5)
polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
        col=alpha("darkgreen", 0.2), lty=0)


# Save fitte curves.
saveRDS(jags.out.constT, "./MCMC_chains/constant_temp.RDS")
saveRDS(jags.out.fluctmin, "./MCMC_chains/fluctmin_temp.RDS")
saveRDS(jags.out.fluctmax, "./MCMC_chains/fluctmax_temp.RDS")