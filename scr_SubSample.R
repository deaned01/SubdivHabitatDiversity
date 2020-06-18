# A script showing the workflow to sub-sample from a stem-mapped forest plot 
# Data are then used to:
# 1. fit shared spp model to alpha diversity;
# 2. calculate predicted and observed shared spp (zeta diversity) among m samples;
# 3. using modelled zeta diversity, calculate for m samples: 
# -- species accumulation curve (and gamma diversity)
# -- number of spp confined to a single patch (single-patch endemics)

# Code to achieve 1-3 above are in script entitled 
load(file="bci_2005.RData")
# contains 3 column dataframe containing data for Barro Colorado Island 50-ha stem-mapped forest plot* 
# data for 2005 census 211845 living stems, 301 species.

# *Available at: https://repository.si.edu/handle/10088/20925; DOI: https://doi.org/10.5479/data.bci.20130603

head(bci_2005)
# spcode= name of species, gx = x coordinate (metres); gy = y coordinate (m)
# plot dimensions = 50 ha, 1000 m in x direction, 500 in y direction.

bci.sad <- as.matrix(unlist(table(bci_2005$spcode)))

# source function to sample from BCI 
source("fn_sampleQuad.R")

args(qr.sample) # function (data, nquad, size, overlap = F, giveup = 20000)
# see function script for explanations
# code to sample 20 quadrats of total area 200 m^2
set.seed(108)
qarea = 200 # change to desired size 
qsize =sqrt(qarea) 
nsamps = 20  # change to desired number of samples (order of zeta)

samp.lst <- list()
reps <- 100
for(i in 1:reps){
  samp.lst[[i]] <- qr.sample(data=bci_2005, nquad=nsamps,size=qsize)
  }

samp20 <- list(samp.lst)
# extract fitting and validation data ####
# mean diversity in 100 qr of qarea

# shared species 
pp1 <- ppfun(samp20, iters=100)

# extract values for calibration ####
shared.spp <- pp1$shared[[1]][,1]
spe <- pp1$spe[[1]]
spp.acc <- pp1$sac[[1]][,1]

# the value to use for fitting c-parameter is mean alpha diversity
alpha.div <- shared.spp[1,1]

# also total spp
gamma.div <- spp.acc[20,1]
# Alternatively, can calculate both directly from the sub-samples
mean(unlist(lapply(samp.lst, function(x) mean(ncol(x))))); gamma.div
mean(unlist(lapply(samp.lst, function(x) apply(x>0,1,sum)))); alpha.div
# should be the same values


args(fitsac)
