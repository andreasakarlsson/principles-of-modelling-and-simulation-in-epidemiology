#Sums of exponentials
rm(list=ls(all=TRUE)) 
#r <- -log(1-1/28)
drugTimes <- seq(0,120,by=5)
drugDoses <- rep(1,length(drugTimes))
drugConcentration <- Vectorize(function(time, r = -log(1-1/28))
  sum(ifelse(time >= drugTimes,drugDoses*exp(-r*(time -drugTimes)),0)))
with(list(t=seq(0,120,by=0.01)),plot(t, drugConcentration(t), type="l"))

#With ODE
rm(list=ls(all=TRUE)) 
require(deSolve)
eventdat <- data.frame(var="y", time= seq(0,120,by=5), value=1, method="add")
parms <- c(rate = -log(27/28))
derivs <- function(t, var, parms)
  list(-parms["rate"]*var["y"])
yini <- c(y = 0)
times <- seq(0, 120, by = 0.01)
out <- lsode(func = derivs, y = yini, times = times, parms = parms, events = list(data = eventdat))
plot(out, main = "", xlab="Time (days)", ylab  = "Drug concentration")

#L1
#With ODE
rm(list=ls(all=TRUE)) 
require(deSolve)
eventdat <- data.frame(var="y", time= seq(0,100,by=1), value=1)
parms <- c(rate = -log(10/5))
derivs <- function(t, var, parms)
  list(-parms["rate"]*var["y"])
yini <- c(y = 0)
#times <- seq(0, 120, by = 0.01)
out <- lsode(func = derivs, y = yini, times = times, parms = parms, events = list(data = eventdat))
plot(out, main = "", xlab="Time (days)", ylab  = "Drug concentration")