require(pomp)

## use pomp with a similar (but not identical) format to deSolve::ode().
## NOTE: func <- function(x,t,parms,delta.t,...) NOT function(t,x,parms,dt,...) !!

psim.pomp <- function(y, times, func, parms=NULL, delta.t=1e-2, ...) {
  names(y) <- paste(names(y),".0",sep="")
  params <- c(parms, y)
  data <- data.frame(
  time=times[-1],
  .report=NA # not currently used  
  )
  
  pomp1 <- pomp(
    data=data,
    times="time",
    t0=times[1],
    rprocess=euler.sim(
      step.fun=func,
      delta.t=delta.t,
      ...
    ),
    rmeasure=function(x,t,params,...) 1 # not currently used 
  )
  ret <- simulate(
    pomp1,
    params=params 
  )
  ret
}

competing.proc.sim <- function (x, t, params, delta.t, ...) {
  tstar <- t+delta.t/2 # mid-point approximation
  rates <- c(params["c1"]*tstar/10, params["c2"]) # intensities
  causes <- reulermultinom(n=1,size=x[1],rate=rates,dt=delta.t)
  `names<-`(x+c(-sum(causes),causes),names(x))
}

psim1 <- psim.pomp(y=c(X=100,Y=0,Z=0),
                   times=seq(0,10,by=0.1),
                   func=competing.proc.sim,
                   parms=c(c1=0.5,c2=0.1)
                   )

with(as(psim1,"data.frame"),     
     matplot(time,cbind(X,Y,Z), type="l", col=1:3))

