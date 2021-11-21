# code to run the SLOSS type simulation (Fig. 5 main text)
load("validationPlots.RData")
source("functions.R")

names(val.plots)
# [1] "bAgg035" "bAgg050" "bRP"     "bReg010" "bci"     "bstrees" "lntrees" "lstrees"

head(val.plots$bAgg035)
# save trees for each validation plot to an object 
btrees <- val.plots$bci
ahtrees <- val.plots$bAgg035
rptrees <- val.plots$bRP
rgtrees <- val.plots$bReg010
bstrees <- val.plots$bstrees
lntrees <- val.plots$lntrees
lstrees <- val.plots$lstrees

# create SAD from each validation plot
bsad <- unlist(table(btrees$spcode))
ahsad <- unlist(table(ahtrees$spcode))
rpsad <- unlist(table(rptrees$spcode))
rgsad <- unlist(table(rgtrees$spcode))
sadbs <- unlist(table(bstrees$spcode))
sadln <-unlist(table(lntrees$spcode))
sadls <-unlist(table(lstrees$spcode))

# bust it up into different amounts of subdiv
# prepare scenarios ####
tota <- 1000*500 # i.e., 50 ha

rem.frac <- seq(from=0.5, to = 0.1, by=-0.1) 
rem.area <- tota*rem.frac

# prepare data holders 
nfrag <- c(1, 2, 4, 8, 16, 32)

mat <- matrix(NA, nrow= length(nfrag), ncol=length(rem.area))
colnames(mat) <- rem.frac
rownames(mat) <- nfrag
for(i in 1:length(nfrag)){
  n <- nfrag[i]
  for(j in 1:length(rem.area)){
    mat[i,j] <- rem.area[j]/n
  }}

# predict the c parameter from the scaling relationships
scens <- data.frame(area = mat[,4])
scens$m <- rownames(mat)
scens$apr <- scens$area/400
scens$cEmp <- 1.06*scens$apr^0.28 # slightly different to the samples used for validation
scens$cAgh <- 0.157*scens$apr^0.778 # 035 fitted equation as per empirical (Eq 8 main text)
scens$cReg <-   seq(from = 1/-1.001, to = 1/-19.999, length.out=6) # regular is approximate b/c power law is not a great fit...

# processing fns ###
it.fn <- function(pr, m) {
  ss <- numeric()
  for(i in 1:m){
    ss[i] <- sum(pr^i)
  }
  return(ss)}

scen.fn <- function(cvec = NULL, avec, mvec, sad, rp = FALSE){
  gam <- spe <- bd <-  numeric(); shsp <- list()
  for(i in 1:length(avec)){
    if(rp)
      {
      pr.calc <- prob.rp(a = avec[i], tota = 500000, sad=bsad)
      } else {
      pr.calc <- prob.fnb(N = sad, cpar = cvec[i], a = avec[i], A = 500000) } 
    
    if(mvec[i] == 1)
      {
      gam[i] <- sum(pr.calc)
      next()
      } else {
        ss <- it.fn(pr.calc, mvec[i])
        spe[i] <- spe.fn(ss)
        bd[i] <- 1-ss[2]/ss[1] 
        shsp[[i]] <- ss
        gam[i] <- gam.fn(ss)
      }
    }
  out <- list(gam = gam, bd = bd, spe=spe, shared=shsp)
  return(out)
  }

# use 10% area ####
emp10 <- scen.fn(cvec = scens$cEmp, avec = scens$area, mvec = scens$m, sad = bsad)

# gamma, spe, bd calcs
e.gd <- (emp10$gam[2:6]-emp10$gam[1])/emp10$gam[1]
e.spe <- emp10$spe[2:6]/emp10$gam[2:6]
e.bd <- emp10$bd[2:6]

# check vals
range(c(emp10$spe[2:6]))
# 28.61215 48.74352
range(c(emp10$gam[2:6]))
# 243.5243 253.0237

# rp 10%
rp10 <- scen.fn(avec = scens$area,mvec = scens$m, sad = bsad, rp=TRUE)
r.gd <- (rp10$gam[2:6]-rp10$gam[1])/rp10$gam[1]
r.bd <- rp10$bd[2:6]
r.spe <- rp10$spe[2:6]/rp10$gam[2:6]

# high agg 10% 
agh10 <- scen.fn(cvec = scens$cAgh, avec = scens$area, mvec = scens$m, sad = ahsad)
a.gd <- (agh10$gam[2:6]-agh10$gam[1])/agh10$gam[1]
a.spe <- agh10$spe[2:6]/agh10$gam[2:6]
a.bd <- agh10$bd[2:6]

# regular pattern 10% 
reg10 <- scen.fn(cvec = scens$cReg, avec = scens$area, mvec = scens$m, sad = bsad)
rg.gd <- (reg10$gam[2:6]-reg10$gam[1])/reg10$gam[1]
rg.bd <- reg10$bd[2:6]
rg.spe <- reg10$spe[2:6]/reg10$gam[2:6]

# relative change calcs 
(emp10$gam[6]-emp10$gam[5])/emp10$gam[1]  # 32 vs 16 patches 
(emp10$gam[3]-emp10$gam[2])/(emp10$gam[6]-emp10$gam[5]) # 4 vs 2 patches
divs <- c(2,4,8,16,32)
sups <-  agh10$gam[2:6]-agh10$gam[1:5] # one subdiv change
# 6.964250 6.111420 4.977772 3.804188 2.730707
plot(divs, sups) 

# package up
library(ggplot2)
library(grid)
mv <- scens$m[2:6]

bd = c(e.bd,r.bd, a.bd, rg.bd)
gamd = c(e.gd, r.gd, a.gd, rg.gd)
spe = c(e.spe, r.spe,a.spe,rg.spe)
sim= c(rep("emp",5),rep("rp",5),rep("agg",5),rep("reg",5))
mv = as.numeric(rep(mv,4))
met <- c(rep('gam',20),rep("spe",20),rep("bd",20))
sims10 <- data.frame(m = rep(mv,3), div = c(gamd,spe,bd), sim=rep(sim,3), metric=met)
sims10$metric <- as.factor(sims10$metric)
sims10$metric <- factor(sims10$metric, levels = c("gam","bd","spe"))
sims10$sim <- as.factor(sims10$sim)
sims10$sim <- factor(sims10$sim, levels= c("emp","rp","agg","reg"),labels=c("Empirical", "Random", "Aggregated", "Regular"))
labs1 <- data.frame(x=c(1,1,1),y=c(0.12,0.3,0.27), lab=c("(a)","(b)","(c)"), metric = as.factor(c("gam","bd","spe")))

# sads ####
empbs <- scen.fn(cvec = scens$cEmp, avec = scens$area, mvec = scens$m, sad = sadbs)
empbs$gam
ebs.gd <- (empbs$gam[2:6]-empbs$gam[1])/empbs$gam[1]
ebs.spe <- empbs$spe[2:6]/empbs$gam[2:6]
ebs.bd <- empbs$bd[2:6]

empln <- scen.fn(cvec = scens$cEmp, avec = scens$area, mvec = scens$m, sad = sadln)
empln$gam
eln.gd <- (empln$gam[2:6]-empln$gam[1])/empln$gam[1]
eln.spe <- empln$spe[2:6]/empln$gam[2:6]
eln.bd <- empln$bd[2:6]

empls <- scen.fn(cvec = scens$cEmp, avec = scens$area, mvec = scens$m, sad = sadls)
empls$gam
els.gd <- (empls$gam[2:6]-empls$gam[1])/empls$gam[1]
els.spe <- empls$spe[2:6]/empls$gam[2:6]
els.bd <- empls$bd[2:6]

mv <- scens$m[2:6]

bd1 = c(e.bd, ebs.bd, eln.bd, els.bd)
gamd1 = c(e.gd, ebs.gd, eln.gd, els.gd)
spe1 = c(e.spe,ebs.spe,eln.spe,els.spe)

sim1= c(rep("emp",5),rep("bs",5),rep("ln",5),rep("ls",5))
mv1 = as.numeric(rep(mv,4))
met1 <- c(rep('gam',20),rep("spe",20),rep("bd",20))
sads10 <- data.frame(m = rep(mv1,3), div = c(gamd1,spe1,bd1), sim=rep(sim1,3), metric=met1)
sads10$metric <- as.factor(sads10$metric)
sads10$metric <- factor(sads10$metric, levels = c("gam","bd","spe"))
sads10$sim <- as.factor(sads10$sim)
sads10$sim <- factor(sads10$sim, levels= c("emp","bs","ln","ls"), labels=c("Empirical", "Broken stick", "Lognormal", "Log series"))
labs2 <- data.frame(x=c(1,1,1),y=c(0.12,0.3,0.27), lab=c("(d)","(e)","(f)"), metric = as.factor(c("gam","bd","spe")))

# combo panel ####
dim(sims10)
dim(sads10)
levels(pldat$sim)
pldat <- rbind(sims10,sads10)
pldat$type <- c(rep("agg", nrow(sims10)),rep("sad",nrow(sads10)))
pldat$type = as.factor(pldat$type)
pldat$type <- factor(pldat$type, levels = c("agg","sad"), labels=c("Vary spatial pattern", "Vary SAD"))
names(pldat) # "m"      "div"    "sim"    "metric" "type"   "leg"   

labs3 <- data.frame(x=c(1,1,1,1,1,1),y=c(0.12,0.3,0.3,0.12,0.3,0.3), lab=c("(a)","(b)","(c)","(d)","(e)","(f)"), 
                    metric = as.factor(c("gam","bd","spe","gam","bd","spe")),
                    type=as.factor(c("agg","agg","agg","sad","sad","sad")))

labs3$type <- factor(labs3$type, levels = c("agg","sad"), labels=c("Vary spatial pattern", "Vary SAD"))
pldat$leg <- as.factor(paste(pldat$type, pldat$sim))

levels(as.factor(pldat$leg))
library(ggplot2)

comfig <- ggplot( data = pldat, aes(x = m, y = div))+
  theme_minimal()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_blank(),
        legend.spacing.x = unit(0.6, "cm"),
        #legend.key.width =unit(4, "cm"),
        legend.position = "bottom",#c(0.84, 0.8),
        legend.title=element_blank(),
        panel.grid = element_blank())+
  
  facet_grid(metric~type, scales="free_y")+
  ylab("")+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab(bquote("Subdivision (number of patches)"))+
  scale_color_manual(values=c("#0072B2", "#009E73","#E69F00", "#999999", "#000000", "#F0E442",  "#D55E00","#CC79A7"), guide=guide_legend(nrow=2, byrow= FALSE))+#, name="Spatial pattern \n   SAD")+
  geom_line(aes(x=m,y=div, col=sim))+
  geom_point(aes(),pch=21, size=2, fill="white")+ 
  geom_text(aes(x=x,y=y,label = lab),fontface="bold",data = labs3); print(comfig)

g2 <- ggplotGrob(comfig)
yax <- which(g2$layout$name=="ylab-l")
g2[["grobs"]][[yax]]$children[[1]]$label <- c('Prop. SPE','Pairwise beta','Change in total species')
g2[["grobs"]][[yax]]$children[[1]]$y <- grid::unit(seq(0.15, 0.85, length=3), "npc")
grid.draw(g2)


