#  PROJECT 
#  =======
#  Principles of Modelling and Simulation in Epidemiology
#  Laboratory exercise 3
#
#  DESCRIPTION
#  ===========
#  Suggested solutions

# AUTHOR AND DATE
# Andreas Karlsson, November 2013

## @knitr Requirements
require(deSolve) #For the deterministic solutions (also for initial and environmental stochasticity)
require(GillespieSSA) #Gillespie Stochastic Simulation Algorithm with "Explicit tau-leap"
## More info at http://www.jstatsoft.org/v25/i12/paper
require(ggplot2) #Extravagant plotting tool (not necessary to answer the questions)
require(plyr) #Handy data manipulation tool (not necessary to answer the questions)
require(parallel) # For Parallel Computation (not necessary to answer the questions)
require(data.table) # Similar to a data.frame but faster (not necessary to answer the questions)

## =============================================================================
## 3 Non-infection models
## =============================================================================
## @knitr Stage_Comp_conf
rm(list=ls(all=TRUE))
mc_ssa <- function(i) {
  init <- c(X = rbinom(n = 1, size = 75, prob =0.15))
  parms <- c(T = 144)
  a <- c("X / T")
  nu <- matrix(c(+1),ncol=1)
  out <- ssa(init, a , nu, parms, tf=53, tau = 0.1, method="ETL", maxWallTime=5, simName='Cancer Progression')
  #'Explicit tau-leap' => user-defined step size
  if (i ==1)
    ssa.plot(out)
  colnames(out$data)[1] <- 'Time'
  Repetition <- rep(i,nrow(out$data))
  Sample <- 1:nrow(out$data)
  rep_out <- data.table(Repetition, Sample, Time = out$data[,1], X = out$data[,2], 
                        Model = rep('Single Compartment',nrow(out$data)))
  return(rep_out)
}
numWorkers <- 3
system.time(res <- mclapply(1:1000, mc_ssa, mc.cores = numWorkers)) # ~4 x faster than a loop
system.time(sim_out <- rbindlist(res)) #Same as do.call("rbind",res), but faster

within_summary <- ddply(sim_out, .(Repetition), summarise, Progression = tail(X, n=1))
percent_prog <- sum(within_summary$Progression!=0)/nrow(within_summary) * 100
print(percent_prog)

rm(list=ls(all=TRUE))
mc_ssa <- function(i) {
  init <- c(SITU_1 = rbinom(n = 1, size = 75, prob =0.15),SITU_2 =0, SITU_3 =0, SITU_4 =0, CANCER_4 =0)
  parms <- c(T1 = 36)
  a <- c("SITU_1 / T1","SITU_2 / T1","SITU_3 / T1","SITU_4 / T1")
  nu <- matrix(c(-1, 0, 0, 0,
                 +1,-1, 0, 0,
                  0,+1,-1, 0,
                  0, 0,+1,-1,
                  0, 0, 0,+1), nrow=5, byrow=T)
  out <- ssa(init, a , nu, parms, tf=53, tau = 0.1, method="ETL", maxWallTime=5, simName='Multiple Compartments')
  #'Explicit tau-leap' => user-defined step size
  if (i ==1)
    ssa.plot(out)
  colnames(out$data)[1] <- 'Time'
  Repetition <- rep(i,nrow(out$data))
  Sample <- 1:nrow(out$data)
  rep_out <- data.table(Repetition, Sample, Time = out$data[,1], X = out$data[,5], 
                        Model = rep('Multiple Compartments',nrow(out$data)))
  return(rep_out)
}
numWorkers <- 3
system.time(res <- mclapply(1:1000, mc_ssa, mc.cores = numWorkers)) # ~4 x faster than a loop
system.time(sim_out <- rbindlist(res)) #Same as do.call("rbind",res), but faster

within_summary <- ddply(sim_out, .(Repetition), summarise, Progression = tail(X, n=1))
percent_prog <- sum(within_summary$Progression!=0)/nrow(within_summary) * 100
print(percent_prog)

## =============================================================================
## 4 Infection models
## =============================================================================
## @knitr det_SIR
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -r * S * I
    dI <- r * S * I - I / T1
    dR <- I / T1    
    return(list(c(dS, dI, dR)))
  })
}
init <- c(S = 1e3, I = 1, R = 0.0)
parameters <- c(r = 1e-4, T1 =12)
times <- seq(0, 1e3, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
out$time <- NULL
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)

rm(list=ls(all=TRUE))
siir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -r * S * ( I1 + I2 )
    dI1 <- - dS - I1 / (T1 / 2)
    dI2 <- (I1 - I2) / (T1 / 2)
    dR <-  I2 / (T1 / 2)    
    return(list(c(dS, dI1, dI2, dR)))
  })
}

init <- c(S = 1e3, I1 = 1, I2 = 0, R = 0)
parameters <- c(r = 1e-4, T1 =12)
times <- seq(0, 1e3, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = siir, parms = parameters))
out$time <- NULL
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIIR Model", lwd = 1, lty = 1, bty = "l", col = 2:5)
legend("topright", c("Susceptibles", "I1", "I2", "Recovereds"), pch = 1, col = 2:5)


rm(list=ls(all=TRUE))
siiir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -r * S * ( I1 + I2 +I3 )
    dI1 <- -dS - I1 / (T1 / 3)
    dI2 <- (I1 - I2) / (T1 / 3)
    dI3 <- (I2 - I3) / (T1 / 3)
    dR <-  I3  / (T1 /3)
    return(list(c(dS, dI1, dI2, dI3, dR)))
  })
}

init <- c(S = 1e3, I1 = 1, I2 = 0, I3 = 0, R = 0)
parameters <- c(r = 1e-4, T1 =12)
times <- seq(0, 1e3, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = siiir, parms = parameters))
out$time <- NULL
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIIIR Model", lwd = 1, lty = 1, bty = "l", col = 2:6)
legend("topright", c("Susceptibles", "I1", "I2", "I3", "Recovereds"), pch = 1, col = 2:6)

rm(list=ls(all=TRUE))
si5r <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -r * S * ( I1 + I2 +I3 + I4 + I5)
    dI1 <- -dS - I1 / (T1 / 5)
    dI2 <- (I1 - I2) / (T1 / 5)
    dI3 <- (I2 - I3) / (T1 / 5)
    dI4 <- (I3 - I4) / (T1 / 5)
    dI5 <- (I4 - I5) / (T1 / 5)
    dR <-  I5 / (T1 / 5)  
    return(list(c(dS, dI1, dI2, dI3, dI4, dI5, dR)))
  })
}

init <- c(S = 1e3, I1 = 1, I2 = 0, I3 = 0, I4 = 0, I5 = 0, R = 0)
parameters <- c(r = 1e-4, T1 =12)
times <- seq(0, 1e3, by = 0.1)
out <- as.data.frame(ode(y = init, times = times, func = si5r, parms = parameters))
out$time <- NULL
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SI5R Model", lwd = 1, lty = 1, bty = "l", col = 2:8)
legend("topright", c("Susceptibles", "I1", "I2", "I3", "I4", "I5", "Recovereds"), pch = 1, col = 2:8)


## @knitr stoc_SIR
rm(list=ls(all=TRUE))
mc_ssa <- function(i) { 
  parms <- c(r=1e-4, T1 = 12)
  init <- c(S = 1000, I = 1, R = 0)
  a <- c("r * S * I", " I / T1")
  nu <- matrix(c(-1, 0, +1, -1, 0, +1), nrow=3, byrow=T)
  out <- ssa(init, a , nu, parms, tf=10e8, tau = 1, method="ETL", ignoreNegativeState=TRUE)
  #'Explicit tau-leap' => user-defined step size
  Repetition <- rep(i,nrow(out$data))
  Sample <- 1:nrow(out$data)
  rep_out <- data.table(Repetition, Sample, Time = out$data[,1], S = out$data[,2], R = out$data[,4], 
                        Model = rep('SIR',nrow(out$data)))
  return(rep_out)
}
numWorkers <- 3
system.time(res <- mclapply(1:1000, mc_ssa, mc.cores = numWorkers)) # ~4 x faster than a loop
system.time(sim_out <- rbindlist(res)) #Same as do.call("rbind",res), but faster

## Plot results
ggplot(sim_out, aes(x=Time, y=R, colour=Repetition)) + geom_line(aes(group = Repetition))

## @knitr stoc_SIIR
rm(list=ls(all=TRUE))
mc_ssa <- function(i) { 
  parms <- c(r=1e-4, T2 = 12 / 2)
  init <- c(S = 1000, I1 = 1, I2 = 0, R = 0)
  a <- c("r * S * (I1 + I2)", " I1 / T2", " I2 / T2" )
  nu <- matrix(c(-1, 0, 0, +1, -1, 0, 0, +1, -1, 0, 0, +1), nrow=4, byrow=T)
  out <- ssa(init, a , nu, parms, tf=10e8, tau = 1, method="ETL", ignoreNegativeState=TRUE)
  #'Explicit tau-leap' => user-defined step size
  Repetition <- rep(i,nrow(out$data))
  Sample <- 1:nrow(out$data)
  rep_out <- data.table(Repetition, Sample, Time = out$data[,1], S = out$data[,2], R = out$data[,5], 
                        Model = rep('SIIR',nrow(out$data)))
  return(rep_out)
}
numWorkers <- 3
system.time(res <- mclapply(1:1000, mc_ssa, mc.cores = numWorkers)) # ~4 x faster than a loop
system.time(sim_out <- rbindlist(res)) #Same as do.call("rbind",res), but faster

## Plot results
ggplot(sim_out, aes(x=Time, y=R, colour=Repetition)) + geom_line(aes(group = Repetition))


## @knitr stoc_SIIIR
rm(list=ls(all=TRUE))
mc_ssa <- function(i) { 
  parms <- c(r=1e-4, T3 = 12 / 3)
  init <- c(S = 1000, I1 = 1, I2 = 0, I3 = 0, R = 0)
  a <- c("r * S * (I1 + I2 + I3)", " I1 / T3", " I2 / T3", " I3 / T3")
  nu <- matrix(c(-1, 0, 0, 0, +1, -1, 0, 0, 0, +1, -1, 0, 0, 0,+1, -1, 0, 0, 0, +1), nrow=5, byrow=T)
  out <- ssa(init, a , nu, parms, tf=10e8, tau = 1, method="ETL", ignoreNegativeState=TRUE)
  #'Explicit tau-leap' => user-defined step size
  Repetition <- rep(i,nrow(out$data))
  Sample <- 1:nrow(out$data)
  rep_out <- data.table(Repetition, Sample, Time = out$data[,1], S = out$data[,2], R = out$data[,6], 
                        Model = rep('SIIIR',nrow(out$data)))
  return(rep_out)
}
numWorkers <- 3
system.time(res <- mclapply(1:1000, mc_ssa, mc.cores = numWorkers)) # ~4 x faster than a loop
system.time(sim_out <- rbindlist(res)) #Same as do.call("rbind",res), but faster

## Plot results
ggplot(sim_out, aes(x=Time, y=R, colour=Repetition)) + geom_line(aes(group = Repetition)) +
  facet_wrap( ~ Model, ncol=2)

## @knitr stoc_SI5R
rm(list=ls(all=TRUE))
mc_ssa <- function(i) {
  parms <- c(r=1e-4, T5 = 12 / 5)
  init <- c(S = 1000, I1 = 1, I2 = 0, I3 = 0, I4 = 0, I5 = 0, R = 0)
  a <- c("r * S * (I1 + I2 + I3 + I4 + I5)", " I1 / T5", " I2 / T5", " I3 / T5"," I4 / T5"," I5 / T5")
  nu <- matrix(c(-1,  0,  0,  0,  0,  0,
                 +1, -1,  0,  0,  0,  0,
                 0, +1, -1,  0,  0,  0,
                 0,  0, +1, -1,  0,  0,
                 0,  0,  0, +1, -1,  0, 
                 0,  0,  0,  0, +1, -1,
                 0,  0,  0,  0,  0, +1), nrow=7, byrow=T)
  out <- ssa(init, a , nu, parms, tf=10e8, tau = 0.1, method="ETL", ignoreNegativeState=TRUE)
  Repetition <- rep(i,nrow(out$data))
  Sample <- 1:nrow(out$data)
  rep_out <- data.table(Repetition, Sample, Time = out$data[,1], S = out$data[,2], R = out$data[,8], 
                        Model = rep('SI5R',nrow(out$data)))
  return(rep_out)
}
numWorkers <- 3
system.time(res <- mclapply(1:1000, mc_ssa, mc.cores = numWorkers)) # ~4 x faster than a loop
system.time(sim_out <- rbindlist(res)) #Same as do.call("rbind",res), but faster

## Plot results
ggplot(sim_out, aes(x=Time, y=R, colour=Repetition)) + geom_line(aes(group = Repetition))
## =============================================================================

## @knitr poisson_system
pos <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dLS <- +r * 1000 - LS / T1
    dLS1 <- +r * 1000 - LS1 / T2
    dLS2 <- (LS1 - LS2) / T2
    dLS3 <- (LS2 - LS3) / T2
    dLSP1 <- +r * 1000 + (- LSP1 - LSP3) / T3
    dLSP2 <- (LSP1 - LSP3) / T3
    dLSP3 <- (LSP1 + LSP2 - LSP3) / T3
    return(list(c( dLS, dLS1, dLS2, dLS3, dLSP1, dLSP2, dLSP3)))
  })
}
init <- c(LS = 0, LS1 = 0, LS2 = 0, LS3 = 0, LSP1 = 0, LSP2 = 0, LSP3 = 0)
parameters <- c(r = 0.0025, T1 = 6, T2 = 2, T3 = 3)
times <- seq(0, 100, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = pos, parms = parameters))
out$time <- NULL
out$LS1 <- NULL
out$LS2 <- NULL
out$LSP1 <- NULL
out$LSP2 <- NULL
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "Poisson Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
matplot(times, cumsum(out)/times, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "Poisson Model", lwd = 1, lty = 1, bty = "l", col = 2:4)

#legend("topright", c("Susceptibles", "I1", "I2", "I3", "Recovereds"), pch = 1, col = 2:6)
















