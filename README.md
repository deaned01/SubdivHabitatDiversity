## SubdivHabitatDiversity
This is a quick introduction to the models presented in “A null model for quantifying the geometric effect of habitat subdivision on species diversity”.  

The R code used in the paper are in file ‘functions.R’ and the aim here is to illustrate their use. The R-code provides a brief explanation of the arguments as well.  

Basically there are two parts to the paper, a validation against empirical data and a theoretical exploration of the implications. I’ll give an example in the reverse order here.

### Estimating diversity
First we need some data - if we assume a value for the community scaling parameter, we only need the species abundance distribution (SAD). 
Here I’ve used data from BCI - that is the Barro Colorado Island forest dynamics plot. I use the 2005 census, with 211845 living stems and 301 species. To get this:  

load('BCI_SAD.RData')

Of course, the SAD could be simulated assuming any number of different abundance distributions…

The models are all in ‘functions.R’.
source("functions.R")
args(fitSS.NB)
function (area, sad, cpar, m, tota = 5e+05) 

The arguments are common to the functions and represent:
–area = sampling grain - in the same units as tota.
–sad = species abundance distribution over the study extent
–cpar = community level scaling parameter at sampling grain ‘area' - this scales roughly as a power function of sampling grain   
–tota = extent over which the number of individuals in the species abundance distribution were counted (defaults to 50 ha)

To predict the diversity of a set of samples from we first calculate zeta diversity and then plug it into the formula derived in Hui & McGeoch (2014).  

Zeta can be calculated using the negative binomial or finite NB versions. This is the neg bin:  
z20 <- fitSS.NB(area=200, sad=bciSAD, cpar = 0.88, m = 20)  

Then once you have zeta, can calculate the total number of species (akin to gamma diversity) using gam.fn()    

gam.fn(z20)  
[1] 141.8386   

Or the number of species in only a single one of those 20 patches using spe.fn()  

spe.fn(z20)  
[1] 45.9334  

There's also bd.fn() for Sorensen dissimilarity.  

### Comparing observed and predicted
1. negative binomial
For comparison with empirical data I’ve included some data obtained from BCI.  

The object ‘calstats’ contains values calculated from 100 repeat samples of 20 randomly positioned 20 x 20 m quadrats.  

load("calstats.RData")
names(calstats)
[1] "zeta"  "spe"   "beta"  "gamma"  

We have four objects, each giving the mean value from the 100 samples from BCI along with 95% sampling intervals:

–‘zeta’ is zeta diversity (average number of species shared) in 1, 2, … , 20 samples.
–‘spe’ is the mean number of species found in only 1 quadrat (single patch endemics)
–‘beta’ is the mean Sorensen dissimilarity
–‘gamma’ is the total species richness

We can predict zeta using the non-random shared species equation and then use this to predict the diversity patterns from the zeta components (this is all derived in Hui & McGeoch 2014).  

First, we need to estimate the community scaling parameter for the 20 x 20 sampling grain.  

alpha.div = calstats$zeta[1,1];alpha.div # mean number of spp shared in 1 sample = alpha diversity
[1] 49.0575  

Once we have this, we can use function fit.c.NB() to estimate the parameter. (NB: If you want to get a scaling relationship for c, repeat the sampling and fitting steps at a few sampling grains, then fit Eq 8 in the main text).  

cpar.fit <- fit.c.NB(obs = alpha.div, area = 400, sad = bciSAD, tota= 500000, low=0, upp=100)  

cpar.fit$cpar [1] 1.038374

This is the value for the community scaling parameter, c, at this sampling grain. Now we just plug it into fitSS.NB() to estimate zeta diversity in the 20 samples…

zeta.est <- fitSS.NB(area = 400, sad = bciSAD, cpar = cpar.fit$cpar, m = 20, tota = 500000); zeta.est  

To compare with the empirical data:  

plot(1:20, calstats$zeta[,1], type="l",ylim=c(0,60),ylab= 'Zeta diversity', xlab="Order of zeta")  
lines(1:20, calstats$zeta[,2], lty= 2)  
lines(1:20, calstats$zeta[,3], lty= 2)  
points(1:20, zeta.est)  
legend(15,50,legend=c("Predicted","Observed"),pch=c(1,NA),lty=c(NA,1))  

Not perfect, but it’s an approximate model. Let’s see how it does for the diversity metrics.  

First total species richness of the samples:

gam.fn(zeta.est) 
[1] 171.0906  

unlist(calstats$gamma)  
mean   hi95  low95 
175.62 188.00 163.00 

Number of species found in a single patch:
spe.fn(zeta.est) 
[1] 46.0355  

unlist(calstats$spe)  
mean  hi95 low95 
43.96 56.00 35.00 

And Sorensen dissimilarity:
bd.fn(zeta.est)  
[1] 0.4755007

unlist(calstats$beta)
  mean      hi95     low95 
0.4936184 0.6155491 0.3616667 

2. Finite negative binomial  
For completeness, here’s the fitting for the FNB model. It takes a long time and does not seem to provide an improvement that warrants this in my view.  

Of course, I’m pretty agricultural in my R skills, so I’m sure others much smarter than I will be able to work out a means to speed this up a whole lot…

To fit c as with the neg bin:  
fnb.c <- fit.c.fnb(obs = alpha.div, area =400, sad = bciSAD, tota= 500000, low=0, upp=10)$cpar
[1] 1.021112  

And predict zeta using the FNB:     
zeta.est.fnb <- fitSS.FNB(sad = bciSAD, cpar = fnb.c, area = 400, tota = 500000, m = 20)  

Compiled by Dave Deane, Nov 2021, d.deane@latrobe.edu.au

Reference  
Hui, C. and McGeoch, M. (2014) Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns American Naturalist 2014 Vol. 184 Issue 5 Pages 684-694  
