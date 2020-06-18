# Appendix C from Modelling shared species in sub-divided habitat 
# draft manuscript by D.Deane, D. Xing, C. Hui, M.McGeoch and F. He
# Code written and maintained by DC Deane deaned01@gmail.com

# This appendix includes a set of functions to:
# fit the non-random shared species model to mean alpha diversity (i.e., estimate parameter c)
# calculate shared species based on the fitted estimate of c
# calculate gamma diversity for m patches of area a
# calculate species confined to a single patch

# Requires R  package Rmpfr to deal with the large values in the gamma fn
# install.packages("Rmpfr")


fitcpar <- function(obs, area, sad, tota= 500000, low=0.00001, upp=20, c.start=1){
  # Estimates parameter c for the finite negative binomial shared species model (Eq 6) 
  # calls css.fnb (objective function required by the optimizer)
  # 
  # Args:
  #    obs= vector of observed alpha diversity at the sampling area/s of interest
  #    area = vector of area/s corresponding to alpha diversity value/s
  #    tota= total study extent (defaults to 50 ha as for BCI; change if not true)
  #    sad = global species abundance distribution (ie for tota) 
  #    low, upp = parameter search range lower and upper bounds 
  #    c.start = starting value for c to be passed to optimizer
  #
  # Returns: 
  #    Either a named list (if only one value estimated) or, dataframe with 2 columns: 
  #    cpar = estimate for c; fit = logical flag indicating convergence (0=fail, 1= converged)
  
  if(length(area)!=length(obs)){
    cat("\n", "Args: obs and area must have equal length",
        "\n")
    
    len = c(length(obs),length(area))
    nom <- c("obs =","area =")
    return(paste(nom, len)) }
  
  fitvec <- estpar <- numeric() # fitvec = logical flag; estpar = etimated cpar  
  for(i in 1:length(obs)){
    
    fit1 = tryCatch({
      optim(fn = css.fnb, obs = obs[i], par = c(cpar = c.start),  area = area[i], tota = tota, 
            sad = sad, method="Brent", lower = low, upper = upp)
    }, warning = function(w) {
      w
    }, error = function(e) {
      e  })
    if(!is(fit1,"error") & length(fit1$par[1])>0) {
      fitvec[i] <- 1
      estpar[i] <- fit1$par[1]
    } else {
      fitvec[i] <- 0
      estpar[i] <- NA}
    print(estpar[i])
  }
  out <- data.frame(cpar= estpar, fit=fitvec)
  return(out)
}

# css.fnb ####
css.fnb <- function(obs, pars, area, tota=500000, sad){
  # calculate alpha div using fnb fit to mean alpha
  # to be used with numerical optimizer (eg optim() ) see- fitcpar()
  # Arguments:
  # a = area, A = total area, sad = species abundance dist; pars = cpar value for optimizer; obs= observed alpha div
  require("Rmpfr")
  cpar = pars[1]
  N = sad
  if(is.list(obs)) obs <- unlist(obs)
  m = length(obs)
  
  alpha <- area/tota
  pr.sar <- numeric()
  
  for(i in 1:length(N)){
    
    pr.sar[i] <- 1 - as.numeric((gamma(as(N[i]*(1+ 1/cpar - alpha/cpar),"mpfr"))*
                                   gamma(as(N[i]/cpar,"mpfr")))/
                                  (gamma(as(N[i]*(1+ 1/cpar),"mpfr"))*
                                     gamma(as(N[i]/cpar*(1- alpha),"mpfr"))))
    
  }
  
  spp <- sum(pr.sar)
  
  dif <- sum(abs(obs - spp))
  
  return(dif)
}

fnbSSfit <- function(sad, cpar, area, tota=500000, m){
  require(Rmpfr)
  # Fn calculates shared species (zeta diversity) in m samples, of area a
  
  # Arguments:
  # sad = vector of abundances (ie SAD over study extent) 
  # area =  subarea of interest
  # tota = total area
  # cpar = c parameter for that grain fit to mean alpha (from fitcpar)
  # m = number of samples of interest for shared count
  
  # Returns:
  # vector of shared species of length m
  
  pr.sar <- numeric()
  alpha = area/tota
  N = sad
  
  for(i in 1:length(N)){
    pr.sar[i] <- 1 - as.numeric((gamma(as(N[i]*(1+ 1/cpar - alpha/cpar),"mpfr"))*
                                   gamma(as(N[i]/cpar,"mpfr")))/
                                  (gamma(as(N[i]*(1+ 1/cpar),"mpfr"))*
                                     gamma(as(N[i]/cpar*(1- alpha),"mpfr"))))
  }
  sar <- (pr.sar)
  ssp <- c()
  
  for(j in 1:m){
    ssp[j] <- sum(sar^j)
  }
  
  out <- list(ssp)
  return(out)
}

fitspacc <- function(obs=NULL, area, sad, tota=500000, cpar, m, plotres = FALSE){#low =1/100, upp = 20){
  
  # Fn to fit and compare a species accumulation curve to observed data 
  # Uses alpha diversity to fit
  
  # Args: 
  # obs =  spp acc curve data 
  # area = sample area in m2 
  # tota = study extent area (assumes 50ha)
  # sad = global sad
  # cpar = fitted c-parameter from fitcpar
  # m = number of samples 
  # plotres= logical flag = 1 for a plot of the observed and fitted
  
  # Returns: 
  # sacm = vector of species accumulated
  
  # Requires:
  # fnbSSfit()
  
  if(plotres & is.null(obs)){
    cat("\n", "Need observed data to plot results",
        "\n")
    }
  
  reps <- m 
  sacm <- numeric()
  obs <- unlist(obs)
  
  ssv <- fnbSSfit(sad = sad, cpar = cpar, area = area, tota = 500000, m = reps)
  modss <- unlist(ssv)
  sacm[1] <- modss[1]
  for(i in 2:reps){
    sse <- modss[1:i]
    len <- length(sse)
    xp <- 1:len
    nx <- (-1)^(xp+1)*sse * choose(len,xp)
    sacm[i] <- (sum(nx))
  }
  if(plotres) plot(1:reps,sacm); points(1:reps, obs,col=2)
  return(sacm)}

spefun <- function(x){
  # Fn to calculate number of single-patch endemic species from expected shared spp 
  
  # Args: 
  # x is a list of MODELLED (not observed) shared species from 1-to-'m' sites 
  # Handles a list at different grains and or reps 
  
  spe <- numeric()
  for(i in 1:length(x)){
    
    # to calculate Sm = number of spp in m patches
    sser <- unlist(x[[i]])
    lenr <- length(sser)
    xpr <- 1:lenr
    
    # to calculate S(m-1) = number in m-1 patches
    xpr.1 <- 1:(lenr-1)
    sser.1 <- sser[-max(xpr)]
    nxr <- (-1)^(xpr+1) * sser * choose(lenr, xpr)
    nxr.1 <- (-1)^(xpr.1+1) * sser.1 * choose((lenr-1), xpr.1)
    gm <- sum(nxr)
    gm.1 <- sum(nxr.1)
    spe[i] <- choose(lenr, 1)*(gm-gm.1)
  }
  return(spe)
}


fnbSPEfit <- function(N, cpar, a, A, m){
  require(Rmpfr)
  # Fn calculates number of species predicted in only a single patch (single-patch endemics) in m patches of area a
  
  # Arguments:
  # N = vector of abundances (ie SAD over study extent) 
  # a =  subarea of interest
  # A = total area
  # cpar = c parameter for that grain fit to mean alpha diversity (see fitcpar)
  # m = number of samples of interest for shared and endemic count
  
  # Returns:
  # Vector of the number of single-patch endemics in m patches of size a
  
  pr.spe <- c()
  alpha = a/A
  
  k <- (1-alpha)*N/cpar
  
  for(i in 1:length(N)){
    pr.spe[i] <- as.numeric((gamma(as(N[i]+ k[i]/(1-alpha) - k[i],"mpfr"))*
                               gamma(as(k[i]/(1-alpha),"mpfr")))/
                              (gamma(as(N[i]+k[i]/(1-alpha),"mpfr"))*
                                 gamma(as(k[i]/(1-alpha) - k[i],"mpfr"))))
  }
  
  esp <- c()
  
  for(j in 1:m){
    esp[j] <- sum(pr.spe^j)
  }
  
  xp <- 1:m
  nx <- (-1)^(xp+1) * esp * choose(m, xp)
  spemod <- sum(nx)
  
  return(spemod)
}
