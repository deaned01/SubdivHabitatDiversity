# R functions from 'A null model for quantifying the geometric effect of habitat subdivision on species diversity'
# 

#***************************************************************************************************************
# Shared spp (zeta) functions ####
# One can use either the negative binomial (quick, recommended) or 
# finite negative binomial (strictly speaking the correct model,for a finite area, but slow and requiring package 
# to deal with the large combinatorial calcs - numerical differences are trivial)

# 1. negative binomial 
predSS.NB <- function(area, sad, cpar, m, tota = 500000){
  # to fit modelled SS with known cpar (par) and SAD
  
  # Arguments:
  #### area = area at which alpha estimated in same units as tota (usually sampling grain in m2)
  #### sad = species abundance distribution (global)
  #### par = c-par for 'area' sampling grain 
  #### m = desired number of shared samples
  #### tota = global area (default assumes 50-ha fdp plot)
  
  # Returns:
  #### vector of zeta diversity of length m 
  
  pr.area = area/tota
  
  S0 <- length(sad)
  kvec <- (sad*pr.area)/cpar
  
  pr.spp <- numeric()
  for(i in 1:S0){
    pr.spp[i] <-  (1 - (1 - pr.area)*(1 + cpar)^-kvec[i])
  }
  
  ssp <- numeric()
  for(j in 1:m){
    ssp[j] <- sum(pr.spp^j)} 
  
  return(ssp)
} 

# 2. Finite negative binomial
# note this can be slow and the improvement in accuracy is unlikely to be worth the wait...
# note too requires package Rmpfr to deal with the large numbers in the gamma function... 
predSS.FNB <- function(sad, cpar, area, tota=500000, m){
  require(Rmpfr)
  # sad = vector of abundances (ie SAD over study extent) 
  # area =  subarea of interest
  # tota = total area
  # cpar = c parameter for that grain fit to mean alpha 
  # m = number of samples of interest for shared and endemic count
  
  # Returns:
  # vector of shared species in m samples
  
  pr.area = area/tota
  
  k <- (pr.area*sad)/cpar
  
  pr.spp <- c()
  for(i in 1:length(sad)){
    pr.spp[i] <- 1 - as.numeric((gamma(as(sad[i]*(1+ 1/cpar - pr.area/cpar),"mpfr"))*
                                   gamma(as(sad[i]/cpar,"mpfr")))/
                                  (gamma(as(sad[i]*(1+ 1/cpar),"mpfr"))*
                                     gamma(as(sad[i]/cpar*(1- pr.area),"mpfr"))))  }
  
  ssp <- c()
  for(j in 1:m){
    ssp[j] <- sum(pr.spp^j)
  }
  
  return(ssp)
}



# 3. random placement

predss.RP <- function(area, tota = 500000, m, sad){
  # a is area of interest
  # A is total study extent for SAD (default value is for 50 ha in metres)
  # returns vector of species shared in 1:m samples
  mvec <- 1:m
  rpmod <- c()
  pr.area <- area/tota
  
  for(i in 1:length(mvec)){
    rpmod[i] <- sum((1-(1- pr.area)^sad)^mvec[i])}
  
  return(rpmod)}


#***************************************************************************************************************
# Diversity functions ####
# Once you have shared species for the relevant setting, the diversity calculations are all done the 
# same way. Three functions are presented here (as used in the paper). 

# All of these functions take one argument: zeta diversity (the mean number of shared spp) in m patches.
# So, the mean number of species in m patches (alpha diversity), the mean number in 2... and so on up to m. 

# See: 
# Hui & McGeoch, Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns American Naturalist 2014 Vol. 184 Issue 5 Pages 684-694

# gamma diversity
gam.fn <- function(x){
  # x is list of shared spp in m patches
  sse <- x
  len <- length(sse)
  xp <- 1:len
  nx <- (-1)^(xp+1) * sse * choose(len, xp)
  gam <- sum(nx)
  return(gam)
}

# beta diversity
# x is list of shared spp in m patches)
# fn returns expected Sorensen pairwise dissimilarity 
bd.fn <- function(x){ 1-x[2]/x[1] } # 

# SPE
# to predict the number of species found in only one patch
spe.fn <- function(x){
  # x is a list of shared spp for m patches
  # fn returns number of spp in only one patch
  sser <- x
  lenr <- length(sser)
  xpr <- 1:lenr
  # to calculate S(m-1)
  xpr.1 <- 1:(lenr-1)
  sser.1 <- sser[-max(xpr)]
  nxr <- (-1)^(xpr+1) * sser * choose(lenr, xpr)
  nxr.1 <- (-1)^(xpr.1+1) * sser.1 * choose((lenr-1), xpr.1)
  gm <- sum(nxr)
  gm.1 <- sum(nxr.1)
  spe <- choose(lenr, 1)*(gm-gm.1)
  return(spe)}


#***************************************************************************************************************
# Validation functions ####
# To test model predictions against actual data, e.g., from a stem-mapped forest
# 1. negative binomial ####

fitc.NB <- function(obs, area, sad, tota= 500000, low=0, upp=100){
  # function to fit the 1-parameter fnb shared species model to a list of observed data
  # basically a wrapper for cparls (see below)
  
  # Args:
  # obs: alpha diversity at grain 'area'
  # area: sampling grain in same units as tota
  # sad: species abundance distribution
  # tota: extent over which sad is calculated
  # low, upp = limits over which to search for c - you will need to adjust to allow negative values 
  # to fit to regular distributions - I find this can be a bit prickly, so several different values might be 
  # necessary to get a fit, try a range like: low = 1/-1.001, upp= 1/-20
  
  #est <- numeric() # fit1 is a logical flag to indicate convergence 
    fit1 = tryCatch({
      optim(fn = cdif.NB, par =  1, obs = obs, area = area, sad = sad, 
            lower = low, upper = upp, method = "Brent")
    }, warning = function(w) {
      w
    }, error = function(e) {
      e  })
    
    if(!is(fit1,"error") & length(fit1$par[1])>0) {
      fit <- 1
      est <- fit1$par[1]
    } else {
      fit <- 0
      est <- NA}
    
  out <- data.frame(cpar= est, fit=fit)
  return(out)
}

cdif.NB <- function(area, sad, par, obs, tota = 500000){
  # helper function to be used with numerical optimizer (eg optim() ) to find the value of 
  # community scaling parameter c that minimises the absolute difference with observed alpha diversity
  
  # Arguments:
  #### area = area at which alpha estimated in same units as tota
  #### sad = species abundance distribution (global)
  #### pars = c-par to fit, 
  #### obs = observed shared species accumulated across the sites represented in the area vector
  #### tota = global area (default assumes 50-ha fdp plot)
  
  # Returns: difference between observed and predicted value for passed c value
  
  pr.area = area/tota
  cpar = as.numeric(par[1])
  
  S0 <- length(sad)
  
  calc <- numeric()
  
  kvec <- sad*pr.area/cpar
  for(i in 1:S0){ 
    calc[i] <-  (1 - (1 - pr.area)*(1 + cpar)^-kvec[i])}
  
  srsim <- sum(calc) 
  dif <- (srsim - obs)
  
  return(abs(dif))
} 

# 2. FNB ####
fitc.FNB <- function(obs, area, sad, tota= 500000, low=-100, upp=100){
  # function to fit the 1-parameter fnb shared species model to a list of observed data
  # basically a wrapper for css.fnb (see below)
  
  # Args:
  # obs: list of observed data. Can be either number of species at a grain or number shared
  # 

  # fit <- est <- numeric() # fit10r is a logical flag to indicate convergence; c20 is estimated cpar at that scale 
 
  fit1 = tryCatch({
      optim(fn = cdif.fnb, obs = obs, par = c(cpar = 1),  area = area, tota = tota, 
            sad = sad, method="Brent", lower = low, upper = upp)
    }, warning = function(w) {
      w
    }, error = function(e) {
      e  })
    
    if(!is(fit1,"error") & length(fit1$par[1])>0) {
      fit <- 1
      est <- fit1$par[1]
    } else {
      fit <- 0
      est <- NA}
  out <- data.frame(cpar= est, fit=fit)
  return(out)
}

cdif.fnb <- function(obs, pars, area, tota = 500000, sad){
  # helper function to be used with numerical optimizer (eg optim() ) to find the value of 
  # community scaling parameter c that minimises the absolute difference with observed alpha diversity

  # Arguments:
  # a = area, tota = total area, sad = species abundance dist.; pars = cpar value for optimizer; obs= observed alpha div
  
  require("Rmpfr")
  cpar = pars[1]
  
  pr.area <- area/tota

  prob.sar <- c()
  
  for(i in 1:length(sad)){
    
    prob.sar[i] <- 1 - as.numeric((gamma(as(sad[i]*(1+ 1/cpar - pr.area/cpar),"mpfr"))*
                                   gamma(as(sad[i]/cpar,"mpfr")))/
                                  (gamma(as(sad[i]*(1+ 1/cpar),"mpfr"))*
                                     gamma(as(sad[i]/cpar*(1- pr.area),"mpfr"))))
    }
  
  ssp <- sum(prob.sar)
  
  dif <- sum(abs(obs - ssp))
  
  return(dif)
}

# species probs (sums gives # spp) 
prob.fnb <- function(N, cpar, a, A = 500000){
  require(Rmpfr)
  # N = vector of abundances (ie SAD over study extent, A) 
  # cpar = c parameter for FNBD model
  # a = area of habitat/sample in same units as A, study extent
  # A = total study extent - default is set to 50 ha
  
  # Returns:
  # vector of probabilities for each species being present in area a 
  #
  # Sum probs to obtain number of species
  # Raise to mth power and sum to obtain number shared in m samples of size a
  
  pr.sar <- c()
  alpha <- a/A
  k <- (alpha*N)/cpar
  for(i in 1:length(N)){
    pr.sar[i] <- 1 - as.numeric((gamma(as(N[i]+k[i]/alpha-k[i],"mpfr"))*
                                   gamma(as(k[i]/alpha,"mpfr")))/
                                  (gamma(as(N[i]+k[i]/alpha,"mpfr"))*gamma(as(k[i]/alpha-k[i],"mpfr"))))}
  
  return(pr.sar)
}

prob.rp <- function(a, tota = 500000, sad){ 
  rpmod <- (1-(1-a/tota)^sad)
  return(rpmod)}

