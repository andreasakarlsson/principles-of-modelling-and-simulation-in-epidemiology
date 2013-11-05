# #Sums of exponentials
# rm(list=ls(all=TRUE)) 
# #r <- -log(1-1/28)
# drugTimes <- seq(0,120,by=5)
# drugDoses <- rep(1,length(drugTimes))
# drugConcentration <- Vectorize(function(time, r = -log(1-1/28))
#   sum(ifelse(time >= drugTimes,drugDoses*exp(-r*(time -drugTimes)),0)))
# with(list(t=seq(0,120,by=0.01)),plot(t, drugConcentration(t), type="l"))
# 
# #With ODE
# rm(list=ls(all=TRUE)) 
# require(deSolve)
# eventdat <- data.frame(var="y", time= seq(0,120,by=5), value=1, method="add")
# parms <- c(rate = -log(27/28))
# derivs <- function(t, var, parms)
#   list(-parms["rate"]*var["y"])
# yini <- c(y = 0)
# times <- seq(0, 120, by = 0.01)
# out <- lsode(func = derivs, y = yini, times = times, parms = parms, events = list(data = eventdat))
# plot(out, main = "", xlab="Time (days)", ylab  = "Drug concentration")
# 
# rm(list=ls(all=TRUE)) 
# require(pomp)
# gompertz.proc.sim <- function (x, t, params, delta.t, ...)
# {
#   ## unpack the parameters:
#   r <- params["r"]
#   K <- params["K"]
#   sigma <- params["sigma"]
#   ## the state at time t:
#   X <- x["X"]
#   ## generate a log-normal random variable:
#   eps <- exp(rnorm(n=1,mean=0,sd=sigma))
#   ## compute the state at time t+delta.t:
#   S <- exp(-r*delta.t)
#   xnew <- c(X=unname(K^(1-S)*X^S*eps))
#   return(xnew)
# }
# 
# 
# gompertz.meas.sim <- function (x, t, params, ...) {
#   ## unpack the parameters:
#   tau <- params["tau"]
#   ## state at time t:
#   X <- x["X"]
#   ## generate a simulated observation:
#   y <- c(Y=unname(rlnorm(n=1,meanlog=log(X),sd=tau)))
#   return(y)
# }
# 
# gompertz.meas.dens <- function (y, x, t, params, log, ...) {
#   ## unpack the parameters:
#   tau <- params["tau"]
#   ## state at time t:
#   X <- x["X"]
#   ## observation at time t:
#   Y <- y["Y"]
#   ## compute the likelihood of Y|X,tau
#   f <- dlnorm(x=Y,meanlog=log(X),sdlog=tau,log=log)
#   return(f)
# }
# 
# gompertz <- pomp(
#   data=data.frame(
#     time=1:100,
#     Y=NA
#   ),
#   times="time",
#   rprocess=discrete.time.sim(
#     step.fun=gompertz.proc.sim,
#     delta.t=1
#   ),
#   rmeasure=gompertz.meas.sim,
#   t0=0
# )
# 
# theta <- c(r=0.1,K=1,sigma=0.1,tau=0.1,X.0=1)
# gompertz <- simulate(gompertz,params=theta)
# plot(gompertz,variables="Y")

#2
#With ODE
rm(list=ls(all=TRUE)) 
require(deSolve)
require(ggplot2)
require(reshape2)

## =============================================================================
## 2.2 Possitive feedback
## =============================================================================
parms <- c(Tc1 = 5, Tc2 = 10, Tc3 = 20)
derivs <- function(t, var, parms)
  list(var[c("T5","T10","T20")]/parms[c("Tc1","Tc2","Tc3")])
yini <- c(T5 = 1, T10=1, T20=1)
times <- seq(0, 10, by = 0.01)
out <- lsode(func = derivs, y = yini, times = times, parms = parms)

##-----------------------------
## Plot results
##-----------------------------
data <- melt(data.frame(out), id="time")
ggplot(data, aes(x=time, y=value, colour=variable)) + geom_line() + 
  labs(title = 'Exponential Increase', y = 'AU')


## =============================================================================
## 2.3a Negative feedback
## =============================================================================
rm(list=ls(all=TRUE)) 
parms <- c(Tc1 = 5, Tc2 = 10, Tc3 = 20)
derivs <- function(t, var, parms)
  list(-var[c("T5","T10","T20")]/parms[c("Tc1","Tc2","Tc3")])
yini <- c(T5 = 100, T10=100, T20=100)
times <- seq(0, 50, by = 0.01)
out <- lsode(func = derivs, y = yini, times = times, parms = parms)

##-----------------------------
## Plot results
##-----------------------------
data <- melt(data.frame(out), id="time")
ggplot(data, aes(x=time, y=value, colour=variable)) + geom_line() + 
  labs(title = 'Exponential Decrease', y = 'AU')

## =============================================================================
## 2.3b Avarage Sojurn Time
## =============================================================================
rm(list=ls(all=TRUE)) 
parms <- c(Tc = 10)
derivs <- function(t, var, parms)
  list(-var["Y"] / parms[c("Tc")])
yini <- c(Y = 100)
times <- seq(0, 100, by = 0.01)
out <- data.frame(lsode(func = derivs, y = yini, times = times, parms = parms))

# Integrating with by summing
Avarage_Time <- 1 / yini * sum(out$Y) * (max(times) / length(times))

##-----------------------------
## Plot results
##-----------------------------
data <- melt(data.frame(out), id="time")
ggplot(data, aes(x=time, y=value, colour=variable)) + geom_line() + 
  ggtitle('Avarage Sojourn Time') +  ylab('AU') +geom_vline(xintercept = Avarage_Time)


## =============================================================================
## 2.4 Negative feedback over two compartments
## =============================================================================
rm(list=ls(all=TRUE)) 
oscillation <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX2 <- -X1
    dX1 <- X2    
    return(list(c(dX1, dX2)))
  })
}

init <- c(X1 = 20, X2 = 0)
parameters <- NULL
times <- seq(0, 10, by = 0.01)
out <- as.data.frame(ode(y = init, times = times, func = oscillation, parms = parameters))
out$time <- NULL

#par(mfrow=c(1,2),pty = "s")
par(mfrow=c(2,2))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
matplot(times, out, type = "l", xlab = "Time", ylab = "AU", main = "Negative feedback over two compartments", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend("topright", c("X1", "X2"), pch = 1, col = 2:4)
plot(out$X1,out$X2,'l', ylab='X2', xlab='X1')



## =============================================================================
## 2.5 Dispertion by a dynamic structure
## =============================================================================

SPCmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    import <- sigimp(t)
    dS <- import - b*S*P + g*C #substrate
    dP <- c*S*P - d*C*P #producer
    dC <- e*P*C - f*C #consumer
    res <- c(dS, dP, dC)
    list(res)
  })
}
## Parameters
parms <- c(b = 0.0, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)
## vector of timesteps
times <- seq(0, 100, length = 101)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                             import = rep(0,length(times))))
signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2
sigimp <- approxfun(signal$times, signal$import, rule = 2)

## Start values for steady state
y <- xstart <- c(S = 1, P = 1, C = 1)
## Solving
out <- lsoda(xstart, times, SPCmod, parms)
## Plotting
mf <- par("mfrow")
plot(out, main = c("substrate", "producer", "consumer"))
plot(out[,"P"], out[,"C"], type = "l", xlab = "producer", ylab = "consumer")
par(mfrow = mf)


## =============================================================================
## 2.5 Dispertion by a dynamic structure
## =============================================================================
rm(list=ls(all=TRUE)) 
deriv <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    import <- sigimp(time)
    dX1 <- (import - X1)/T0
    dX2 <- (X1 - X2)/T0
    dX3 <- (X2 - X3)/T0 
    return(list(c(dX1, dX2,dX3)))
  })
}

init <- c(X1 = 0, X2 = 0, X3 = 0)
parameters <- c(T0=10)
dt <- 0.01
times <- seq(0, 100, by = dt)

## external signal with triangle impulse with area 100
signal <- as.data.frame(list(times = times, import = rep(0,length(times))))
signal$import[signal$times >= 0 & signal$times <= 0] <- 100 / dt * 2 #The 2 since it is a triangle
sigimp <- approxfun(signal$times, signal$import, rule = 2)

out <- as.data.frame(ode(y = init, times = times, func = deriv, parms = parameters))
out$time <- NULL

##-----------------------------
## Plot results
##-----------------------------
par(mfrow=c(1,1))
matplot(times, cbind(sigimp(times),out), type = 'l', xlab = 'Time', ylab = 'AU', main = 'Dispersion of a flow that passes a dynamic structure', lwd = 1, lty = 1, bty = 'l', col = 2:5, ylim=c(0,12))
legend('topright', c('Pulse','X1', 'X2','X3'), lty = 1, col = 2:5)


## =============================================================================
## 2.6 A logistic and an SI structure
## =============================================================================
library(deSolve)
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

init <- c(S = 1-1e-6, I = 1e-6, 0.0)
parameters <- c(beta = 1.4247, gamma = 0.14286)
times <- seq(0, 70, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
out$time <- NULL

##-----------------------------
## Plot results
##-----------------------------
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)

amat <- diag(4)
diag(amat) <- NA
svar.a <- SVAR(var.2c, estmethod = "logLik", Amat = amat)
irf(svar.a, impulse = "e", response = c("prod", "rw", "U"), boot =
      FALSE)

require(vars)
data(Canada)
## For VAR
var.2c <- VAR(Canada, p = 2, type = "const")
irf.1 <- irf(var.2c, impulse = "e", response = c("prod", "rw", "U"), boot = FALSE)
plot(irf.1)

require(pomp)
require(deSolve)
sir.proc.sim <- function (x, t, params, delta.t, ...) {
  ## unpack the parameters
  N <- params["N"] # population size
  gamma <- params["gamma"] # recovery rate
  mu <- params["mu"] # birth rate = death rate
  beta <- params["beta"] # contact rate
  foi <- beta*x["I"]/N # the force of infection
  trans <- c(
    rpois(n=1,lambda=mu*N*delta.t), # births are assumed to be Poisson
    reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t), # exits from S
    reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t), # exits from I
    reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t) # exits from R
  )
  ## now connect the compartments
  x[c("S","I","R","cases")]+c(
    trans[1]-trans[2]-trans[3],
    trans[2]-trans[4]-trans[5],
    trans[4]-trans[6],
    trans[4] # accumulate the recoveries
  )
}
simulate(
  pomp(
    data=data.frame(
      time=seq(1/52,15,by=1/52),
      reports=NA
    ),
    times="time",
    t0=0,
    rprocess=euler.sim(
      step.fun=sir.proc.sim,
      delta.t=1/52/20
    ),
    measurement.model=reports~binom(size=cases,prob=rho
    ),
    initializer=function(params, t0, ic.pars, comp.names, ...){
      x0 <- c(S=0,I=0,R=0,cases=0)
      N <- params["N"]
      fracs <- params[ic.pars]
      x0[comp.names] <- round(N*fracs/sum(fracs))
      x0
    },
    zeronames=c("cases"), # 'cases' is an accumulator variable
    ic.pars=c("S0","I0","R0"), # initial condition parameters
    comp.names=c("S","I","R") # names of the compartments
  ),
  params=c(
    N=50000,
    beta=60,gamma=8,mu=1/50,
    rho=0.6,
    S0=8/60,I0=0.002,R0=1-8/60-0.001
  ),
  seed=677573454L
) -> sir
plot(sir)




## =============================================================================
## A simple delay differential equation  
## dy(t) = -y(t-1) ; y(t<0)=1 
## =============================================================================

##-----------------------------
## the derivative function
##-----------------------------
derivs <- function(t, y, parms) {
  if (t < 1)
    dy <- -1
  else
    dy <- - lagvalue(t - 1)
  list(c(dy))
}

##-----------------------------
## initial values and times
##-----------------------------
yinit <- 1
times <- seq(0, 30, 0.1)

##-----------------------------
## solve the model  
##-----------------------------
yout <- dede(y = yinit, times = times, func = derivs, parms = NULL)

##-----------------------------
## display, plot results
##-----------------------------
plot(yout, type = "l", lwd = 2, main = "dy/dt = -y(t-1)")

## =============================================================================
## The infectuous disease model of Hairer; two lags.
## example 4 from Shampine and Thompson, 2000
## solving delay differential equations with dde23
## =============================================================================

##-----------------------------
## the derivative function
##-----------------------------
derivs <- function(t,y,parms) {
  if (t < 1)
    lag1 <- 0.1
  else 
    lag1 <- lagvalue(t - 1,2)
  if (t < 10)
    lag10 <- 0.1
  else 
    lag10 <- lagvalue(t - 10,2)
  
  dy1 <- -y[1] * lag1 + lag10
  dy2 <-  y[1] * lag1 - y[2]
  dy3 <-  y[2] - lag10
  list(c(dy1, dy2, dy3))
}

##-----------------------------
## initial values and times
##-----------------------------
yinit <- c(5, 0.1, 1)
times <- seq(0, 40, by = 0.1)

##-----------------------------
## solve the model  
##-----------------------------
system.time(
  yout <- dede(y = yinit, times = times, func = derivs, parms = NULL)
)

##-----------------------------
## display, plot results
##-----------------------------
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, lty = 1,
        main = "Infectuous disease - Hairer")

## =============================================================================
## time lags + EVENTS triggered by a root function
## The two-wheeled suitcase model 
## example 8 from Shampine and Thompson, 2000
## solving delay differential equations with dde23
## =============================================================================

##-----------------------------
## the derivative function
##-----------------------------
derivs <- function(t, y, parms) {
  if (t < tau)
    lag <- 0
  else 
    lag <- lagvalue(t - tau)
  
  dy1 <- y[2]
  dy2 <- -sign(y[1]) * gam * cos(y[1]) +
    sin(y[1]) - bet * lag[1] + A * sin(omega * t + mu)
  list(c(dy1, dy2))
}

## root and event function
root <- function(t,y,parms) ifelse(t>0, return(y), return(1))
event <- function(t,y,parms) return(c(y[1], y[2]*0.931))

gam = 0.248; bet = 1; tau = 0.1; A = 0.75
omega = 1.37; mu = asin(gam/A)

##-----------------------------
## initial values and times
##-----------------------------
yinit <- c(y = 0, dy = 0)
times <- seq(0, 12, len = 1000)

##-----------------------------
## solve the model  
##-----------------------------
## Note: use a solver that supports both root finding and events, 
##       e.g. lsodar, lsode, lsoda, adams, bdf
yout <- dede(y = yinit, times = times, func = derivs, parms = NULL,  
             method = "lsodar", rootfun = root, events = list(func = event, root = TRUE))

##-----------------------------
## display, plot results
##-----------------------------

plot(yout, which = 1, type = "l", lwd = 2, main = "suitcase model", mfrow = c(1,2))
plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)