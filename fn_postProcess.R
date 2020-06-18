ppfun <- function(dat, iters){
  # extracts validation data from set of samples (or list of sets)
  # dat = output from fn qr.sample
  
  sac <- shared <- sads <- spe <- list() # lists at level of output (elements = length(dat))
  # sads <- avegam <- sr <- avess <- avesac <- avesr <- gamhi <- gamlo <- list()
  
  for(i in 1:length(dat)){
    x <- dat[[i]]
    m <- nrow(x[[1]])
    #browser()
    # declare holders for each iteration 
    sac.obs <- ss.obs <- matrix(NA, nrow = length(x), ncol=m)
    srj <- spej <- gamj <-  numeric()  
    sadj <- list()
    for(j in 1:length(x)){
      
      y <- x[[j]]
      
      zed <- which(apply(y, 2, sum)==0) # check for zero abundance species and remove
      if(length(zed>0)) y <- y[,-zed]
      
      gamj[j] <- ncol(y) # gamma diversity 
      
      sadj[[j]] <- apply(y,2,sum)
      
      ofd <- apply(y>0,2,sum) # occupancy frequency distribution
      spej[j] = length(which(ofd==1)) # number of single patch endemics (found in only one sample)
      
      srj[j] <- mean(apply(y>0,1,sum))
      
      # observed shared species and SAC estimated from 'iters' iterations 
      ssiter <- matrix(NA, nrow = iters, ncol=m)
      saciter <- matrix(NA, nrow =iters, ncol=m)
      
      for(p in 1:iters){
        y <- y[sample(1:nrow(y)),]
        saciter[p,1:(m-1)] <- rowSums(apply(y[1:(m-1), ], 2, cumsum) > 0)
        for(k in 2:m){
          #browser()
          ssiter[p,k] <- length(which(apply(y[1:k,] > 0, 2, sum) == k))
        }
        ssiter[p,1] <- srj[j]
        saciter[p,m] <- gamj[j]
        }
      
      sac.obs[j,] <- apply(saciter, 2, mean)
      ss.obs[j,] <- apply(ssiter,2,mean)
    }
    #browser()
    ssave <- apply(ss.obs, 2, mean)
    sshi <- apply(ss.obs, 2, function(x) quantile(x, 0.975))
    sslo <- apply(ss.obs, 2, function(x) quantile(x, 0.025))
    sacave <- apply(sac.obs, 2, mean)
    sachi <- apply(sac.obs, 2, function(x) quantile(x, 0.975))
    saclo <- apply(sac.obs, 2, function(x) quantile(x, 0.025))
    sac[[i]] <- data.frame(avesac=sacave, hici = sachi, loci = saclo)
    shared[[i]] <-  data.frame(ssave = ssave, hici = sshi, loci = sslo)
    sads[[i]] <- sadj
    avespe <- mean(spej)
    spehi <- quantile(spej, 0.98)
    spelo <- quantile(spej, 0.03)
    spe[[i]] <- data.frame(speave = avespe, hici=spehi, loci=spelo)
    print(i)
  }
  
  out <- list(shared = shared, sac=sac, spe=spe, sads=sads)
  return(out)
}
