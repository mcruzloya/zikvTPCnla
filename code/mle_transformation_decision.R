#library('R2jags')
#library('mcmcplots')
#library('MCMCvis')
library('scales')

setwd("/Users/cruzloya/git/zikvTPCnla/")

data.constT <- read.csv("./data/processed_data/constantT.csv")
data.fluctmax <- read.csv("./data/processed_data/fluctmax.csv")
data.fluctmin <- read.csv("./data/processed_data/fluctmin.csv")

View(data.constT) 
#View(data.constT.old)

#order(data.constT)
#data.constT <- data.constT[with(data.constT, order(data.constT$temp, data.constT$replicate)), ]
#data.constT.old <- data.constT.old[with(data.constT.old, order(data.constT.old$temp,
                                                          data.constT.old$replicate)), ]

#data.constT.old$vir.rate <- data.constT.old$vir.rate / 10^9
#View(data.constT)
#View(data.constT.old)

#data.constT$vir.rate - data.constT.old$vir.rate

data.fluctmin$vir.rate
data.constT$vir.rate

set.seed(42)

par(mfrow=c(1,3))
plot(data.constT$temp, data.constT$vir.rate, xlab="temperature",
     ylab="viruses / hour", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 12*10^9), pch=20, col="steelblue")
plot(data.fluctmin$temp, data.fluctmin$vir.rate, xlab="temperature", 
     ylab="viruses / hour", 
     main="2.5C fluctuation", xlim=c(15, 38), ylim=c(0, 12*10^9), pch=20, col="steelblue")
plot(data.fluctmax$temp, data.fluctmax$vir.rate, xlab="temperature", 
     ylab="viruses / hour", 
     main="5C fluctuation", xlim=c(15, 38), ylim=c(0, 12*10^9), pch=20, col="steelblue")

par(mfrow=c(1,3))
k <- 7
plot(data.constT$temp, log10(1 + data.constT$vir.rate / 10^k), xlab="temperature",
     ylab="log(1 + viral titer / hour)", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")
plot(data.fluctmin$temp, log10(1 + data.fluctmin$vir.rate / 10^k), xlab="temperature", 
     ylab="log(1 + viral titer / hour", 
     main="2.5C fluctuation", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")
plot(data.fluctmax$temp, log10(1 + data.fluctmax$vir.rate / 10^k), xlab="temperature", 
     ylab="log(1 + viral titer / hour)", 
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


nll <- function(theta, ydata, Tdata) {
  Tmin <- theta[1]
  Tmax <- theta[2]
  rmax <- theta[3]
  alpha <- theta[4]
  beta <- theta[5]
  sigma2 <- theta[6]
  k <- theta[7]
  hi <- Inf
  
  # Check parameters in bounds
  if(rmax < 0) {
    return(hi)
  } else if((alpha < 0) | (alpha > 1)) {
    return(hi)
  } else if((beta < 0) | (beta > 1)) {
    return(hi)
  } else if(sigma2 < 0) {
    return(hi)
  } else if((k < 0) | (k > 9)) {
    return(hi)
  }
  
  y_k <- log10(1 + ydata / 10^k)
  mu <- flexTPC(Tdata, Tmin, Tmax, rmax, alpha, beta)
  result <- sum((y_k - mu)^2 / sigma2 + log(10^k + ydata) ) + length(ydata) * log(sqrt(2*pi*sigma2))
  return(result)
}

par(mfrow=c(1,1))
k <- 7
plot(data.constT$temp, log10(1 + data.constT$vir.rate / 10^k), xlab="temperature",
     ylab="log(1 + viral titer / hour)", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")
temps <- seq(15, 40, 0.01)
theta0 <- c(15, 35, 2.5, 0.7, 0.3, 0.25, 7)
lines(temps, flexTPC(temps, theta0[1], theta0[2], theta0[3], theta0[4], theta0[5]))

nll(theta0, ydata=data.constT$vir.rate,
        Tdata=data.constT$temp)

data.constT$vir.rate
data.constT$temp

theta_hat <- optim(theta0, nll, ydata=data.constT$vir.rate,
                   Tdata=data.constT$temp, hessian=TRUE)
theta_hat

theta_hat$par
nll(theta_hat$par, ydata=data.constT$vir.rate,
        Tdata=data.constT$temp)
k
k <- theta_hat$par[7]
k

plot(data.constT$temp, log10(1 + data.constT$vir.rate / 10^k), xlab="temperature",
     ylab="log(1 + viral titer / hour)", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 5), pch=20, col="steelblue")
temps <- seq(15, 40, 0.01)
lines(temps, flexTPC(temps, theta_hat$par[1], theta_hat$par[2],
                     theta_hat$par[3], theta_hat$par[4], theta_hat$par[5]))



