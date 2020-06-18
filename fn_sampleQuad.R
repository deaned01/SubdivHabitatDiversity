# Suite of functions to sample from a 50-ha plot 
# Adapted from code written by Fangliang He, University of Alberta

qr.sample <- function (data, nquad, size, overlap = F, giveup = 20000)   {
  # function to sample from a stem-mapped forest plot, or similar data, giving spp xy location and identity. 
  # Sampling assumes square quadrats of a number and size specified by the user.
  
  # Arguments: 
  #  data = a long-form data frame with x and y coordinates names "gx" and "gy" and species column called "spcode"
  #  nquad = desired number of quadrats to sample 
  #  size = the size of the quadrats (one side of a square qr)
  #  overlap = logical argument specifying whether quadrats can be overlapping 
  #  giveup = threshold number of attempts to position non-overlapping quadrats within sampling space 
  # (function terminates with error message if unable to find solution)
  # 
  # Returns:
  # a matrix of dimensions nquad rows x number of species sampled
  
  # specify dimensions
  x = data$gx
  y = data$gy
  sp <- data$spcode
  
  minx = min(x)
  maxx = max(x)
  miny = min(y)
  maxy = max(y)
  
  # call helper function to position quadrats and plot
  quadrat.xy = quad.fn(minx, maxx, miny, maxy, nquad, size, overlap, giveup)
  xlo = quadrat.xy$x - size/2
  xup = quadrat.xy$x + size/2
  ylo = quadrat.xy$y - size/2
  yup = quadrat.xy$y + size/2
  
  trees <- data.frame(x = data$gx, y = data$gy, spcode=data$spcode)
  
  # specify holder matrix dimensions
  nspp <- length(unique(trees$spcode))
  sppnams <- levels(trees$spcode)
  
  z1 <- spp1 <- list() # list to hold individual 
  
  for(j in 1:nquad){
    sppI <- which(trees$x >= xlo[j] & trees$x <= xup[j] & 
                    trees$y >= ylo[j] & trees$y <= yup[j])
    
    spp1[[j]] <- trees$spcode[sppI]
    z1[[j]] <- as.numeric(unlist(table(spp1[[j]])))
    }
  sppDat = do.call(rbind, z1)
  colnames(sppDat) <- sppnams
  zed <-  which(apply(sppDat,2,sum)==0) 
  # remove columns for species with zero abundance
  if(length(zed)!=0) sppDat <- sppDat[,-zed]
  return(sppDat)
  } # End function 


# quadrat count main function quadrat.count.main ####
quadrat.count.main <- function (data, nquad, size, overlap = F, giveup = 200) 
  {
  npt = numeric()
  x = data$gx
  y = data$gy
  minx = min(x)
  maxx = max(x)
  miny = min(y)
  maxy = max(y)
  quadrat.xy = quad.fn(minx, maxx, miny, maxy, nquad, size, overlap, giveup)
  xlo = quadrat.xy$x - size/2
  xup = quadrat.xy$x + size/2
  ylo = quadrat.xy$y - size/2
  yup = quadrat.xy$y + size/2
  
  for (i in 1:nquad) {
    npt[i] = length(x[(x >= xlo[i] & x < xup[i]) & (y >= ylo[i] & y < yup[i])])
    }
  return(npt)
  }

# quadgen.fn ####
quadgen.fn <- function (size, boundary) 
  {
  quadrat = list(x = runif(1, min = boundary[1] + size/2, max = boundary[2] - 
                             size/2), y = runif(1, min = boundary[3] + size/2, max = boundary[4] - size/2))
  return(quadrat)
}

# quadencr.fn ####
quadencr.fn <- function (element, tlist, size) 
  {
  answer = F
  for (i in 1:length(tlist$x)) {
    dist = list(x = abs(element$x - tlist$x[i]), y = abs(element$y - tlist$y[i]))
    if (max(dist$x, dist$y) <= size) 
      answer = T
    }
  answer
  }

# quad.fn #####
quad.fn <- function (minx, maxx, miny, maxy, nquad, size, overlap = F, giveup = 200) 
  {
  boundary = c(minx, maxx, miny, maxy)
  if(nquad==1){
    quadlist = quadgen.fn(size, boundary)
    } else {
      quadlist = quaddadd.fn(nquad, size, overlap, giveup, boundary)
    }
  return(quadlist)
}

# quaddadd.fn #### 
quaddadd.fn <- function (nquad, size, overlap, giveup, boundary, addlist) 
  {
  if (missing(addlist)) {
    addlist = quadgen.fn(size, boundary)
    needone = 0
    } else {
      needone = 1
      }
  for (i in 2:(nquad + needone)) {
    gen = quadgen.fn(size, boundary)
    count = 0
    if (overlap == F) {
      while (quadencr.fn(gen, addlist, size) == T) {
        gen = quadgen.fn(size, boundary)
        count = count + 1
        if (count > giveup) {
          cat("\n", "Search for non-overlapping quadrats reached the giveup level (program terminated)",
             "\n")
          return(list(x = 0, y = 0))
        }
      }
    }
    addlist <- list(x = cbind(addlist$x, gen$x), y = cbind(addlist$y,gen$y))
    }
  return(data.frame(x = array(addlist$x), y = array(addlist$y)))
}
