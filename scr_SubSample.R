# A script showing the workflow to sub-sample from a stem-mapped forest plot 
# Data are then used to:
# 1. fit shared spp model to alpha diversity;
# 2. calculate predicted and observed shared spp (zeta diversity) among m samples;
# 3. using modelled zeta diversity, calculate for m samples: 
# -- species accumulation curve (and gamma diversity)
# -- number of spp confined to a single patch (single-patch endemics)

# Code to achieve 1-3 above are in script entitled 
load(file="bci_2005.RData")
# contains 3 column dataframe containing data for Barro Colorado Island 50-ha stem-mapped forest plot
# Available at: https://repository.si.edu/handle/10088/20925
# DOI: https://doi.org/10.5479/data.bci.20130603

head(bci_2005)
# spcode= name of species, gx = x coordinate (metres); gy = y coordinate (m)
# plot dimensions = 50 ha, 1000 m in x direction, 500 in y direction.


