## 
# This script processes the raw data to calculate viral replication rates.

setwd("/Users/cruzloya/git/zikvTPCnla/")

## Read experiment with multiple rates
data <- list()

temps <- c(18, 20, 22, 28, 30, 32)
dt_vec <- vector(mode="integer")
constant <- vector(mode="numeric")
fluctmin <- vector(mode="numeric")
fluctmax <- vector(mode="numeric")

for(i in 1:length(temps)) {
  temperature <- temps[i]
  filename <- paste("./data/raw_data/tpcdata_", temperature, ".csv", sep="")
  curr_data <- read.csv(filename)
  data[[i]] <- curr_data
  
  # Find time interval from this experiment.
  nt <- length(curr_data$time_days)
  dt <- curr_data$time_days[nt] - curr_data$time_days[nt-1]
  dt_vec <- c(dt_vec, dt)
}

## Plot some examples from raw data.

par(mfrow=c(1,3))
Tidx = 1
plot(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_1, pch=20,
     xlab="time [hours]", ylab="viral titer", main="18+-2.5°C, r1")
plot(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_2, pch=20,
     xlab="time [hours]", ylab="viral titer", main="18+-2.5°C, r2")
plot(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_3, pch=20,
     xlab="time [hours]", ylab="viral titer", main="18+-2.5°C, r3")


colors = c("steelblue", "firebrick", "forestgreen")

Tidx = 1 
plot(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_1, pch=20,
     xlab="time [days]", ylab="viral titer", main="18+-2.5°C",
     col=colors[1])
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_1, col=colors[1])
points(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_2, col=colors[2],
       pch=20)
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_2, col=colors[2])
points(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_3, col=colors[3],
       pch=20)
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_3, col=colors[3])

Tidx = 3 
plot(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_1, pch=20,
     xlab="time [days]", ylab="viral titer", main="22+-2.5°C",
     col=colors[1])
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_1, col=colors[1])
points(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_2, col=colors[2],
       pch=20)
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_2, col=colors[2])
points(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_3, col=colors[3],
       pch=20)
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmin_3, col=colors[3])

Tidx = 4
plot(data[[Tidx]]$time_days, data[[Tidx]]$fluctmax_1, pch=20,
     xlab="time [days]", ylab="viral titer", main="28+-5°C", col=colors[1],
     ylim=c(0, max(c(data[[Tidx]]$fluctmax_1, data[[Tidx]]$fluctmax_2, data[[Tidx]]$fluctmax_3))))
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmax_1, col=colors[1])
points(data[[Tidx]]$time_days, data[[Tidx]]$fluctmax_2, col=colors[2],
       pch=20)
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmax_2, col=colors[2])
points(data[[Tidx]]$time_days, data[[Tidx]]$fluctmax_3, col=colors[3],
       pch=20)
lines(data[[Tidx]]$time_days, data[[Tidx]]$fluctmax_3, col=colors[3])





processed_rates <- data.frame(temp=rep(rep(temps, 3), 3), dt=rep(rep(dt_vec, 3), 3),
                              replicate=rep(c(rep(1, 6), rep(2, 6), rep(3, 6)), 3),
                              condition=(c(rep("baseline", 18), rep("fluctmin", 18),
                                          rep("fluctmax", 18))),
                              temp.fluct=(c(rep(0, 18), rep(2.5, 18),
                                            rep(5, 18)))
                              )


processed_rates$peak.time <- NaN
processed_rates$peak.titer <- NaN

processed_rates
#processed_rates['condition']
data[[1]]
for(i in 1:6) {
  temp <- temps[i]
  curr_data <- data[[i]]
  times <- curr_data$time_days
  for(condition in c("baseline", "fluctmin", "fluctmax")) {
    for(replicate in 1:3) {
      titers <- curr_data[paste(condition, replicate, sep="_")]
      peak.titer <- max(titers)
      peak.timeidx <- which.max(as.numeric(unlist(titers)))
      peak.time <- times[peak.timeidx]
      print(paste(condition, temp, replicate, peak.time) )
      processed_rates$peak.titer[(processed_rates["condition"] == condition) &
                      (processed_rates["replicate"] == replicate) &
                      (processed_rates["temp"] == temp)  ] <- peak.titer
      processed_rates$peak.time[(processed_rates["condition"] == condition) &
                        (processed_rates["replicate"] == replicate) &
                        (processed_rates["temp"] == temp)] <- peak.time
    }
  }
}
processed_rates
processed_rates$vir.rate <- 24 * (processed_rates$peak.titer / processed_rates$dt)
processed_rates$inf.rate <- 1 / processed_rates$peak.time



processed_rates
data[[1]]
head(curr_data, 10)
#View(curr_data)

data.constT <- subset(processed_rates, processed_rates$condition == "baseline")
data.fluctmin <- subset(processed_rates, processed_rates$condition == "fluctmin")
data.fluctmax <- subset(processed_rates, processed_rates$condition == "fluctmax")


data.constT
## Add data from newest experiment with constant temperatures
data.newexp.raw <- read.csv("./data/raw_data/additional_constant_temps.csv")

  

data.newexp <- data.frame(temp=c(rep(24, 3), rep(26, 3)), dt=2,
                          replicate=rep(1:3, 2),
                          condition=rep("baseline", 6),
                          temp.fluct=rep(0, 6))
data.constT
data.newexp$peak.time <- NaN
data.newexp$peak.titer <- NaN
data.newexp$vir.rate <- NaN
data.newexp$inf.rate <- NaN

data.newexp
temps <- c(24, 26)
n.temps <- 2
for(i in 1:n.temps) {
  curr_temp <- temps[i]
  for(rep in 1:3) {
    print(curr_temp)
    print(rep)
    print(data.newexp.raw$temp == curr_temp)
    curr_data <- subset(data.newexp.raw, ((data.newexp.raw$temp == curr_temp) &
                          (data.newexp.raw$replicate == rep)))
    print(curr_data)
    titers <- curr_data$virus
    peak.titer <- max(titers)
    peak.timeidx <- which.max(as.numeric(unlist(titers)))
    peak.time <- curr_data$time_days[peak.timeidx]
    data.newexp$peak.titer[(data.newexp$temp == curr_temp) & (data.newexp$replicate == rep)] <- peak.titer
    data.newexp$peak.time[(data.newexp$temp == curr_temp) & (data.newexp$replicate == rep)] <- peak.time
  }
  
  
  
}

data.newexp
data.newexp$vir.rate <- 24 * (data.newexp$peak.titer /  data.newexp$dt)
data.newexp$inf.rate <- 1 / data.newexp$peak.time

data.newexp


#data.newexp <- read.csv("constantT_processed_newexp.csv")

data.newexp
data.constT

rbind(data.constT, data.newexp)
data.constT <- rbind(data.constT, data.newexp)

data.constT

#View(processed_rates)
par(mfrow=c(1,3))
plot(data.constT$temp, data.constT$vir.rate / 10^9, xlab="temperature",
     ylab="viruses / hour * 10^-9", 
     main="constant temperature", xlim=c(15, 38), ylim=c(0, 12), pch=20, col="steelblue")


plot(data.fluctmin$temp, data.fluctmin$vir.rate / 10^9, xlab="temperature", 
     ylab="viruses / hour * 10^-9", 
     main="2.5C fluctuation", xlim=c(15, 38), ylim=c(0, 12), pch=20, col="steelblue")
plot(data.fluctmax$temp, data.fluctmax$vir.rate / 10^9, xlab="temperature", 
     ylab="viruses / hour * 10^-9", 
     main="5C fluctuation", xlim=c(15, 38), ylim=c(0, 12), pch=20, col="steelblue")


par(mfrow=c(1,3))
plot(data.constT$temp, data.constT$inf.rate, xlab="temperature",
     ylab="1 / peak time", 
     main="constant temperature", xlim=c(15, 38),  ylim=c(0, 0.5), pch=20, col="steelblue")

plot(data.fluctmin$temp, data.fluctmin$inf.rate, xlab="temperature", 
     ylab="1 / peak time", 
     main="2.5C fluctuation", xlim=c(15, 38),  ylim=c(0, 0.5), pch=20, col="steelblue")
plot(data.fluctmax$temp, data.fluctmax$inf.rate, xlab="temperature", 
     ylab="1 / peak time", 
     main="5C fluctuation", xlim=c(15, 38),  ylim=c(0, 0.5), pch=20, col="steelblue")


data.constT
write.csv(data.constT, file="./data/processed_data/constantT.csv")
write.csv(data.fluctmax, file="./data/processed_data/fluctmax.csv")
write.csv(data.fluctmin, file="./data/processed_data/fluctmin.csv")
