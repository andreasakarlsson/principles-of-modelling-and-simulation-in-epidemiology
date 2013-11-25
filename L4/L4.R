#  PROJECT 
#  =======
#  Principles of Modelling and Simulation in Epidemiology
#  Laboratory exercise 4
#
#  DESCRIPTION
#  ===========
#  Suggested solutions

# AUTHOR AND DATE
# Andreas Karlsson, November 2013

## @knitr Requirements
require(deSolve) #For the deterministic solutions (also for initial and environmental stochasticity)
require(stats) #The stats package has the 'optim'-function


## @knitr Defining_input_data
rm(list=ls(all=TRUE))

## Setting the relative tolerance for the optimisations
reltol <- 1e-4

CA_STAT_0_amp = c(0, 0, 0, 0, 1, 4.5, 14, 32, 53, 54, 48, 43, 36.5, 29, 26.5, 24, 22, 18)
CA_STAT_0 <- approxfun(seq(0,85,by=5),CA_STAT_0_amp, rule = 1)

CA_STAT_1_amp <- rep(0,85)
CA_STAT_1_amp[32:55] = c(21.7, 25.4, 29.2, 34.6, 38.7, 42.0, 44.7, 49.4, 54.1, 52.7, 52.3,
              46.9, 41.7, 39.0, 36.8, 32.1, 26.6, 24.5, 22.1, 19.4, 18.1, 20.0,
              19.5, 19.5)
CA_STAT_1 <- approxfun(1:85,CA_STAT_1_amp, rule = 2)

SCR_STAT_amp <- rep(0,85)
SCR_STAT_amp[32:55]  <- c(10.0, 12.8, 16.2, 24.7, 37.5, 54.7, 80.0, 132.0, 178.6, 206.7,
             212.2, 211.4, 202.3, 183.4, 159.9, 123.7, 95.0, 78.0, 61.4, 46.6,
             37.1, 31.3, 29.9, 26.6)
SCR_STAT <- approxfun(1:85,SCR_STAT_amp, rule = 2)

Fraction_alive_amp <- c(1, 0.993, 0.993, 0.992, 0.991, 0.989, 0.987, 0.984, 0.980, 0.973, 0.963, 0.948, 0.925, 0.889, 0.834, 0.746, 0.601, 0.417)
Fraction_alive <- approxfun(seq(0,85,by=5),Fraction_alive_amp, rule = 1)

## @knitr constructing_the_model
CervNoScreen <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    D <- T1 / (3 * P^(1/3) )
    R <- T1 / (3 * (1 - P^(1/3) ))
    dSITU1_0 <- New_insitu(time) - SITU1_0 / R - SITU1_0 / D
    dSITU2_0 <- SITU1_0 / D - SITU2_0 / R - SITU2_0 / D
    dSITU3_0 <- SITU2_0 / D - SITU3_0 / R - SITU3_0 / D
    dInvasive_0 <-  SITU3_0 / D - Invasive_0 / T2
    Ca_diag_0 <- (Invasive_0 / T2 )
    return(list(c(dSITU1_0, dSITU2_0, dSITU3_0, dInvasive_0), Ca_diag_0 = Ca_diag_0))
  })
}
CervScreen <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    D <- T1 / (3 * P^(1/3) )
    R <- T1 / (3 * (1 - P^(1/3) ))
    SITU_1 <- SITU1_1 + SITU2_1 + SITU3_1
    dSITU1_1 <- New_insitu(time) - SITU1_1 / R - SITU1_1 / D - SCR_STAT(time) * (SITU1_1/max(SITU_1, 1e-15))
    dSITU2_1 <- SITU1_1 / D - SITU2_1 / R - SITU2_1 / D - SCR_STAT(time) * (SITU2_1/max(SITU_1, 1e-15))
    dSITU3_1 <- SITU2_1 / D - SITU3_1 / R - SITU3_1 / D - SCR_STAT(time) * (SITU3_1/max(SITU_1, 1e-15))
    dInvasive_1 <-  SITU3_1 / D - Invasive_1 / T2 - Invasive_1 * SCR_STAT(time) / max(SITU_1, 1e-15)
    Ca_diag_1 <- Invasive_1 / T2 
    Ca_diag_scr <- Invasive_1 * SCR_STAT(time) / max(SITU_1, 1e-15)
    Ca_diag_both <- Ca_diag_1 + Ca_diag_scr
    return(list(c(dSITU1_1, dSITU2_1, dSITU3_1, dInvasive_1), Ca_diag_both = Ca_diag_both))
  })
}
callCervModels <- function(parameters, NewInsituLookUpTable = FALSE, make_plots = FALSE){
  if (NewInsituLookUpTable){
    New_insitu <<- approxfun(seq(0, 85, by = 2.5), parameters[grep("Y_tri", names(parameters))], rule = 2)
  }
  else{
    New_insitu <<- approxfun(c(0, parameters['a'], parameters['b'], parameters['c'], 85), c(0, 0, parameters['y'], 0, 0), rule = 2)
  }
  dt <- 1
  times <- seq(0, 85, by = dt)
  
  Ca_diag_0 <- 0
  init_0 <- c(SITU1_0 = 0, SITU2_0 = 0, SITU3_0 = 0, Invasive_0 = 0)
  out_nscr <- as.data.frame(ode(y = init_0, times = times, func = CervNoScreen, parms = parameters))
  out_nscr$time <- NULL
  
  Ca_diag_1 <- 0
  init_1 <- c(SITU1_1 = 0, SITU2_1 = 0, SITU3_1 = 0, Invasive_1 = 0)
  out_scr <- as.data.frame(ode(y = init_1, times = times, func = CervScreen, parms = parameters))
  out_scr$time <- NULL
  
  V_0 <- sum((out_nscr['Ca_diag_0'] - CA_STAT_0(times))^2)
  V_1 <- sum((out_scr[times >= 31.5 & times <= 55.5, 5] - CA_STAT_1(times[times >= 31.5 & times <= 55.5]))^2)
  V <- V_0 + V_1 + ifelse(parameters['P'] <= 0 , 1e7, 0) 
  if (make_plots){
    layout(matrix(1:3, 1, 3, byrow = TRUE))
    matplot(times, cbind(New_insitu(times), SCR_STAT(times)), type = "l", xlab = "Time", ylab = "Subjects", main = "Input Function", lwd = 1, lty = 1, bty = "l", col = 2:4)
    legend('topright', c('New_insitu', 'SCR_STAT'), lty = 1, col = 2:3)
    matplot(times, cbind(CA_STAT_0(times),out_nscr['Ca_diag_0']), type = "l", xlab = "Time", ylab = "Subjects", main = "CervNoScreen", lwd = 1, lty = 1, bty = "l", col = 2:3)
    legend('topright', c('CA_STAT_0', 'Ca_diag_0'), lty = 1, col = 2:3)
    matplot(times, cbind(CA_STAT_0(times),CA_STAT_1(times),out_scr['Ca_diag_both']), type = "l", xlab = "Time", ylab = "Subjects", main = "CervScreen", lwd = 1, lty = 1, bty = "l", col = 2:4)
    legend('topright', c('CA_STAT_0', 'CA_STAT_1', 'Ca_diag_0'), lty = 1, col = 2:4)  
    print('Obj. functions:')
    print(round(data.frame(V_0 = V_0, V_1 = V_1, V = V, row.names = NULL)))
  }
  return(V)
}
## @knitr initial_value_results
parameters <- c(P = 0.1, T1 = 10, T2 = 5, a = 20, b = 33, c = 80, y = 550)
res_value <- callCervModels(parameters, make_plots = TRUE)

## @knitr optimised_results
opt_res <- optim(parameters,  callCervModels, method = "BFGS", control=list( reltol)) 
tmp <- callCervModels(parameters = opt_res$par, make_plots = TRUE)
print(opt_res$par)
print(paste('The obj. function was reduced by:', round((res_value - opt_res$value)/res_value *100),'%'))

## @knitr rep_opt_res_with_lookup_table
triang <- approxfun(c(0, opt_res$par['a'], opt_res$par['b'], opt_res$par['c'], 85), c(0, 0, opt_res$par['y'], 0, 0), rule = 2)
triang_time <- seq(0, 85, by = 2.5)
parameters_lookup <- c(opt_res$par[c('P', 'T1', 'T2')], Y_tri = triang(triang_time))
res_triang <- callCervModels(parameters_lookup, NewInsituLookUpTable = TRUE, make_plots = TRUE)

## @knitr opt_res_with_opt_triangle_Nelder_Mead
opt_tri_res_Simplex <- optim(parameters_lookup,  callCervModels, method = c("Nelder-Mead"), control=list( reltol), NewInsituLookUpTable = TRUE) 
tmp <- callCervModels(parameters = opt_tri_res_Simplex$par, NewInsituLookUpTable = TRUE, make_plots = TRUE)
print(opt_tri_res_Simplex$par)

## @knitr opt_res_with_opt_triangle_BFGS
opt_tri_res_BFGS <- optim(parameters_lookup,  callCervModels, method = c("BFGS"), control=list(reltol), NewInsituLookUpTable = TRUE) 
tmp <- callCervModels(parameters = opt_tri_res_BFGS$par, NewInsituLookUpTable = TRUE, make_plots = TRUE)
print(opt_tri_res_BFGS$par)


## @knitr opt_of_screening
AppliedCervScreen <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    D <- T1 / (3 * P^(1/3) )
    R <- T1 / (3 * (1 - P^(1/3) ))
    SITU <- SITU1 + SITU2 + SITU3
    dSITU1 <- New_insitu(time) - SITU1 / R - SITU1 / D
    dSITU2 <- SITU1 / D - SITU2 / R - SITU2 / D
    dSITU3 <- SITU2 / D - SITU3 / R - SITU3 / D
    dInvasive <-  SITU3 / D - Invasive / T2
    Ca_diag_both <- Invasive / T2
    dSUMCANCER <- Ca_diag_both * Fraction_alive(time)
    return(list(c(dSITU1, dSITU2, dSITU3, dInvasive, dSUMCANCER), Ca_diag_both = Ca_diag_both, SITU = SITU))
  })
}
callAppliedCervModels <- function(ages, opt_parameters, make_plots = FALSE){
  New_insitu <<- approxfun(seq(0, 85, by = 2.5), opt_parameters[grep("Y_tri", names(opt_parameters))], rule = 2)
  dt <- 1
  times <- seq(0, 85, by = dt)
  
  eventdat <- data.frame(var = c("SITU1", "SITU2", "SITU3","Invasive"),
                       time = c(rep(ages['Age1'],4), rep(ages['Age2'],4), rep(ages['Age3'],4)) ,
                       value = rep(0.25,12),
                       method = rep("mult",12))
  
  Ca_diag <- 0
  init <- c(SITU1 = 0, SITU2 = 0, SITU3 = 0, Invasive = 0, SUMCANCER = 0)
  out_appl_nscr <- as.data.frame(ode(y = init, times = times, func = AppliedCervScreen, parms = opt_parameters))
  out_appl_scr <- as.data.frame(ode(y = init, times = times, func = AppliedCervScreen, parms = opt_parameters, events = list(data = eventdat)))  
  if (make_plots){
    layout(matrix(1:2, 1, 2, byrow = TRUE))
    matplot(times, cbind(out_appl_nscr$SITU, out_appl_scr$SITU), type = "l", xlab = "Time", ylab = "Ratios", main = "SITU_0 vs SITU_1", lwd = 1, lty = 1, bty = "l", col = 2:4)
    legend('topright', c('SITU_0', 'SITU_1'), lty = 1, col = 2:3)
    matplot(times, cbind(out_appl_nscr['SUMCANCER'],out_appl_scr['SUMCANCER']), type = "l", xlab = "Time", ylab = "Subjects", main = "Cumulated Cancers", lwd = 1, lty = 1, bty = "l", col = 2:3)
    legend('topleft', c('SUMCANCER_0', 'SUMCANCER_1'), lty = 1, col = 2:3)
  }
  
  REDUCTION <- max(out_appl_nscr['SUMCANCER'], na.rm=T) - max(out_appl_scr['SUMCANCER'], na.rm=T)
  return(REDUCTION)
}
# Choosing the optimised parameters from previously
opt_parameters <- opt_tri_res_BFGS$par

## @knitr noscreening
ages <- c(Age1 = 0, Age2 = 0, Age3 = 0)
tmp <- callAppliedCervModels(ages =ages, opt_parameters, make_plots=TRUE)

## @knitr onescreening
ages <- c(Age1 = 40, Age2 = 0, Age3 = 0)
tmp <- callAppliedCervModels(ages =ages, opt_parameters, make_plots=TRUE)

## @knitr optimised_screening_ages
lower = 0
upper = 85
ages <- c(Age1 = 35, Age2 = 45, Age3 = 55)
opt_ages_L_BFGS_B <- optim(ages,  callAppliedCervModels, method = c("L-BFGS-B"), lower = lower, upper = upper, control=list(reltol, fnscale = -1), opt_parameters = opt_parameters) 
tmp <- callAppliedCervModels(ages = round(opt_ages_L_BFGS_B$par), opt_parameters, make_plots = TRUE)
print(opt_ages_L_BFGS_B$par)