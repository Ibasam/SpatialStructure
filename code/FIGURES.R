#########------------------------- GLOBAL DATA AND PARAMETERS -------------------------#############

nSIMUL=50 #100
npop=15
nYears=40
nInit=10
pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
dat<-dat[-1,]



## SCENARIOS
EXPE<-c(10311111, 30311111, 50311111, 70311111, #0-30% No diversity - simple configuration
        10331111, 30331111, 50331111, 70331111, #0-30% Gradual diversity - simple configuration
        10331211, 30331211, 50331211, 70331211, #0-30% Random diversity - simple configuration
        30331221, 50331221, 70331221,           #10-30% Random diversity - complex configuration (distance)
        30331212, 50331212, 70331212,           #10-30% Random diversity - complex configuration (carrying capacity)
        30331222, 50331222, 70331222)           #10-30% Random diversity - complex configuration (both)

scn_nodiv <- c(1:4)
scn_grad <- c(5:8) 
scn_rand <- c(9:12)
scn_spat_10 <- c(10,13,16,19)
scn_all_div <- c(1:12)
scn_all_spat <- c(10:21)


color_scn <- c("grey80","grey60","grey40","grey20",
               "darkseagreen1","mediumseagreen","forestgreen","darkolivegreen",
               "gold", "darkorange","darkorange3", "lightsalmon4",
               "red","red", "red", 
               "dodgerblue","dodgerblue","dodgerblue",
               "darkviolet", "darkviolet", "darkviolet")

#########--------------------------- PACKAGES AND FUNCTIONS ---------------------------#############
library(visNetwork)
library(reshape2)
library(imputeTS)

#########--------------------------- LOAD RESULTS FROM SIMULATIONS ---------------------------#############

res_1SW_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/PHENOGENOTYPE",EXPE[iEXPE],"_50.RData"))
  #load(paste0("results/PHENOGENOTYPE/PHENOGENOTYPE",EXPE[iEXPE],"_MSW_50.RData"))
  res_1SW_scn[[iEXPE]]<-res_1SW
  #res_MSW_scn[[iEXPE]]<-res_MSW
}

nReturns_scn <- list()
Mig_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/DEMOGRAPHY",EXPE[iEXPE],"_50.RData"))
  nReturns_scn[[iEXPE]]<-nReturns
  Mig_scn[[iEXPE]]<-Mig
}

#optimum
optimum<-rep(0,15)

#########---------------------------- CODE FOR FIGURES ----------------------------#############


#COMPUTE METRICS OF LOCAL TRAIT MISMATCH AND EVOLUTIONARY RATES
ltm_vit<-array(,dim=c(npop,length(EXPE)))
ltm_ratio<-array(,dim=c(npop,length(EXPE)))
ltm_end<-array(,dim=c(npop,length(EXPE)))
ltm_init<-array(,dim=c(npop,length(EXPE)))
cv_init<-array(,dim=c(npop,length(EXPE)))
cv_end<-array(,dim=c(npop,length(EXPE)))

for (pop in 1:15) { 
  geno_1SW_tot<-list()
  for (scn in c(1:length(EXPE))){
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  for (scn in 1:length(EXPE)) {
    median_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    cv_gene<-aggregate(exp(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop]), by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), cv)[1:50,2]
    lw<-loess(median_gene~c(1:50),span=0.2)
    lw2<-loess(cv_gene~c(1:50),span=0.2)
    
    ltm_vit[pop,scn]<- min(which(abs(predict(lw))[10:50] <= (mean(abs(predict(lw))[45:50], na.rm=T)+mean(abs(predict(lw))[45:50], na.rm=T)*0.01)))
    ltm_init[pop,scn]<-mean(abs(predict(lw))[5:10], na.rm=T)
    ltm_end[pop,scn]<-mean(abs(predict(lw))[45:50], na.rm=T)
    
    ltm_ratio[pop,scn]<-abs(mean(predict(lw)[5:10], na.rm=T)-mean(predict(lw)[45:50], na.rm=T))
    ltm_ratio[pop,scn]<-ltm_ratio[pop,scn]/ltm_vit[pop,scn]
    
    cv_init[pop,scn]<-mean(abs(predict(lw2))[5:10], na.rm=T)
    cv_end[pop,scn]<-mean(abs(predict(lw2))[45:50], na.rm=T)
    
  }
}


##############################################################
### Fig.3 BOXPLOT ALL SCN DIVERSITY (SIMPLE CONFIGURATION) ###
##############################################################

##Local trait mismatch
gap<-c(-.2,0,.2)
plot(NULL, xlim=c(0.5,5.5), ylim=c(0,0.14), xlab="Dispersal rate",ylab="Local trait mismatch", xaxt='n')
mtext(c("Init.", "0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:5)
ltm_init_median_nodiv <- apply(ltm_init[,scn_nodiv],1,mean)
ltm_init_median_grad <- apply(ltm_init[,scn_grad],1,mean)
ltm_init_median_rand <- apply(ltm_init[,scn_rand],1,mean)
abline(v=1.5,lty=3)
points(c(1:5+gap[1]),c(median(ltm_init_median_nodiv),apply(ltm_end[,1:4],2,median)),pch=19, cex=1.3)
segments(c(1:5+gap[1]),c(quantile(ltm_init_median_nodiv,probs=.025),apply(ltm_end[,1:4],2,quantile,probs=.025)),
         c(1:5+gap[1]),c(quantile(ltm_init_median_nodiv,probs=.975),apply(ltm_end[,1:4],2,quantile,probs=.975)),lwd=1)
segments(c(1:5+gap[1]),c(quantile(ltm_init_median_nodiv,probs=.25),apply(ltm_end[,1:4],2,quantile,probs=.25)),
         c(1:5+gap[1]),c(quantile(ltm_init_median_nodiv,probs=.75),apply(ltm_end[,1:4],2,quantile,probs=.75)),lwd=2)

points(c(1:5+gap[2]),c(median(ltm_init_median_grad),apply(ltm_end[,5:8],2,median)),pch=19, col="forestgreen", cex=1.3)
segments(c(1:5+gap[2]),c(quantile(ltm_init_median_grad,probs=.025),apply(ltm_end[,5:8],2,quantile,probs=.025)),
         c(1:5+gap[2]),c(quantile(ltm_init_median_grad,probs=.975),apply(ltm_end[,5:8],2,quantile,probs=.975)),lwd=1, col="forestgreen")
segments(c(1:5+gap[2]),c(quantile(ltm_init_median_grad,probs=.25),apply(ltm_end[,5:8],2,quantile,probs=.25)),
         c(1:5+gap[2]),c(quantile(ltm_init_median_grad,probs=.75),apply(ltm_end[,5:8],2,quantile,probs=.75)),lwd=2, col="forestgreen")

points(c(1:5+gap[3]),c(median(ltm_init_median_rand),apply(ltm_end[,9:12],2,median)),pch=19, col="orange", cex=1.3)
segments(c(1:5+gap[3]),c(quantile(ltm_init_median_rand,probs=.025),apply(ltm_end[,9:12],2,quantile,probs=.025)),
         c(1:5+gap[3]),c(quantile(ltm_init_median_rand,probs=.975),apply(ltm_end[,9:12],2,quantile,probs=.975)),lwd=1, col="orange")
segments(c(1:5+gap[3]),c(quantile(ltm_init_median_rand,probs=.25),apply(ltm_end[,9:12],2,quantile,probs=.25)),
         c(1:5+gap[3]),c(quantile(ltm_init_median_rand,probs=.75),apply(ltm_end[,9:12],2,quantile,probs=.75)),lwd=2, col="orange")
legend("topright", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=c("black","forestgreen","orange"), bty="n", ncol=1)


##Evolutionary rate
plot(NULL, xlim=c(0.5,4.5), ylim=c(0,14), xlab="Dispersal rate",ylab="Evolutionary rate", xaxt='n')
mtext(c("0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4)
points(c(1:4+gap[1]),apply(ltm_ratio[,1:4],2,median)*1000,pch=19, cex=1.3)
segments(c(1:4+gap[1]),apply(ltm_ratio[,1:4],2,quantile,probs=.025)*1000,
         c(1:4+gap[1]),apply(ltm_ratio[,1:4],2,quantile,probs=.975)*1000,lwd=1)
segments(c(1:4+gap[1]),apply(ltm_ratio[,1:4],2,quantile,probs=.25)*1000,
         c(1:4+gap[1]),apply(ltm_ratio[,1:4],2,quantile,probs=.75)*1000,lwd=2)

points(c(1:4+gap[2]),apply(ltm_ratio[-8,5:8],2,median)*1000,pch=19, col="forestgreen", cex=1.3)
segments(c(1:4+gap[2]),apply(ltm_ratio[-8,5:8],2,quantile,probs=.025)*1000,
         c(1:4+gap[2]),apply(ltm_ratio[-8,5:8],2,quantile,probs=.975)*1000,lwd=1, col="forestgreen")
segments(c(1:4+gap[2]),apply(ltm_ratio[-8,5:8],2,quantile,probs=.25)*1000,
         c(1:4+gap[2]),apply(ltm_ratio[-8,5:8],2,quantile,probs=.75)*1000,lwd=2, col="forestgreen")

points(c(1:4+gap[3]),apply(ltm_ratio[-6,9:12],2,median)*1000,pch=19, col="orange", cex=1.3)
segments(c(1:4+gap[3]),apply(ltm_ratio[-6,9:12],2,quantile,probs=.025)*1000,
         c(1:4+gap[3]),apply(ltm_ratio[-6,9:12],2,quantile,probs=.975)*1000,lwd=1, col="orange")
segments(c(1:4+gap[3]),apply(ltm_ratio[-6,9:12],2,quantile,probs=.25)*1000,
         c(1:4+gap[3]),apply(ltm_ratio[-6,9:12],2,quantile,probs=.75)*1000,lwd=2, col="orange")

legend("topleft", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=c("black","forestgreen","orange"), bty="n", ncol=1)


################################################################################
### Fig.4 TEMPORAL EVOLUTION OF gG SOME POPULATIONS SCN SIMPLE CONFIGURATION ###
################################################################################

#compute proportion of Immigrants
P.im<-array(,dim=c((nYears+nInit),npop,nSIMUL,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      for (t in 1:(nYears+nInit)) {
        P.im[t,pop,simul,scn] <- Mig_scn[[scn]][[1]][[simul]]$NIm[t,pop] / nReturns_scn[[scn]][[1]][[simul]][t,pop]
      }
    }
  }
}

P.im.median<-array(,dim=c((nYears+nInit),npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (pop in 1:npop) {
    for (t in 1:(nYears+nInit)) {
      P.im.median[t,pop,scn] <- median (P.im[t,pop,,scn], na.rm=T)
    }
  }
}


#graph
par(mfrow=c(3,3))
for (pop in 1:15) {
  geno_1SW_tot<-list()
    for (scn in scn_nodiv){ #define scenarios #scn_grad #scn_rand
      geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
      for (simul in 2:nSIMUL) {
        geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
      }
      geno_1SW_tot[[scn]] <- geno_1SW
   }
  
  color_transparent <- function(color) adjustcolor(color, alpha.f = 0.3)
  gap <- c(-.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3,-.2)
  
  plot(NULL, xlim=c(0,50), ylim=c(-.15,.15), xlab="Years",ylab="Growth potential", main=pop, cex.lab=1.5,cex.axis=1.2)
  rect(0,-.15,10,.15, density = NULL, col = color_transparent("grey"), border = NA)
  abline(h=0, lty=2)
   
  for (scn in scn_nodiv) { #define scenarios here too
      median_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
      lw<-loess(median_gene~c(1:50),span=0.2)
      lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
  
      #if(scn %in% c(2:4,6:8,10:12) ) {
      if(scn >1 ) { #plot immigrants (for scn with dispersal rates > 0)
          median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
          lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
          lines((12:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[10:50,pop,scn]),lty=3)
      }
  }
  legend("topright", cex=.8, legend=c("No div. - No dispersal","No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Gradual div. - No dispersal","Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[5:8], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - No dispersal","Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[9:12], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(1,4,7,10)], bty="n", ncol=1)
  legend("bottomright", cex=.8, legend=c("Philo","Immi"), border=NA,lty=c(1,3), bty="n")
}



###################################################
### Fig.5 BOXPLOT ALL SCN SPATIAL CONFIGURATION ###
###################################################

##Local trait mismatch
gap<-c(-.3,-.1,.1,.3)
plot(NULL, xlim=c(0.5,4.5), ylim=c(0,0.14), xlab="Dispersal rate",ylab="Local trait mismatch", xaxt='n')
mtext(c("Init.","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4)
ltm_init_median_rand <- apply(ltm_init[,c(10:12)],1,mean)
ltm_init_median_dist <- apply(ltm_init[,c(13:15)],1,mean)
ltm_init_median_size <- apply(ltm_init[,c(16:18)],1,mean)
ltm_init_median_both <- apply(ltm_init[,c(19:21)],1,mean)
abline(v=1.5,lty=3)
points(c(1:4+gap[1]),c(median(ltm_init_median_rand),apply(ltm_end[,10:12],2,median)),pch=19, col="orange", cex=1.3)
segments(c(1:4+gap[1]),c(quantile(ltm_init_median_rand,probs=.025),apply(ltm_end[,10:12],2,quantile,probs=.025)),
         c(1:4+gap[1]),c(quantile(ltm_init_median_rand,probs=.975),apply(ltm_end[,10:12],2,quantile,probs=.975)),lwd=1, col="orange")
segments(c(1:4+gap[1]),c(quantile(ltm_init_median_rand,probs=.25),apply(ltm_end[,10:12],2,quantile,probs=.25)),
         c(1:4+gap[1]),c(quantile(ltm_init_median_rand,probs=.75),apply(ltm_end[,10:12],2,quantile,probs=.75)),lwd=2, col="orange")

points(c(1:4+gap[2]),c(median(ltm_init_median_dist),apply(ltm_end[,13:15],2,median)),pch=19, col="red", cex=1.3)
segments(c(1:4+gap[2]),c(quantile(ltm_init_median_dist,probs=.025),apply(ltm_end[,13:15],2,quantile,probs=.025)),
         c(1:4+gap[2]),c(quantile(ltm_init_median_dist,probs=.975),apply(ltm_end[,13:15],2,quantile,probs=.975)),lwd=1, col="red")
segments(c(1:4+gap[2]),c(quantile(ltm_init_median_dist,probs=.25),apply(ltm_end[,13:15],2,quantile,probs=.25)),
         c(1:4+gap[2]),c(quantile(ltm_init_median_dist,probs=.75),apply(ltm_end[,13:15],2,quantile,probs=.75)),lwd=2, col="red")

points(c(1:4+gap[3]),c(median(ltm_init_median_size),apply(ltm_end[,16:18],2,median)),pch=19, col="dodgerblue", cex=1.3)
segments(c(1:4+gap[3]),c(quantile(ltm_init_median_size,probs=.025),apply(ltm_end[,16:18],2,quantile,probs=.025)),
         c(1:4+gap[3]),c(quantile(ltm_init_median_size,probs=.975),apply(ltm_end[,16:18],2,quantile,probs=.975)),lwd=1, col="dodgerblue")
segments(c(1:4+gap[3]),c(quantile(ltm_init_median_size,probs=.25),apply(ltm_end[,16:18],2,quantile,probs=.25)),
         c(1:4+gap[3]),c(quantile(ltm_init_median_size,probs=.75),apply(ltm_end[,16:18],2,quantile,probs=.75)),lwd=2, col="dodgerblue")

points(c(1:4+gap[4]),c(median(ltm_init_median_both),apply(ltm_end[,19:21],2,median)),pch=19, col="darkviolet", cex=1.3)
segments(c(1:4+gap[4]),c(quantile(ltm_init_median_both,probs=.025),apply(ltm_end[,19:21],2,quantile,probs=.025)),
         c(1:4+gap[4]),c(quantile(ltm_init_median_both,probs=.975),apply(ltm_end[,19:21],2,quantile,probs=.975)),lwd=1, col="darkviolet")
segments(c(1:4+gap[4]),c(quantile(ltm_init_median_both,probs=.25),apply(ltm_end[,19:21],2,quantile,probs=.25)),
         c(1:4+gap[4]),c(quantile(ltm_init_median_both,probs=.75),apply(ltm_end[,19:21],2,quantile,probs=.75)),lwd=2, col="darkviolet")

legend("topright", cex=.8, legend=c("Random div.","Random + Dist","Random + Size","Random +D+S"), col=2:3, border=NA,fill=c("orange","red","dodgerblue","darkviolet"), bty="n", ncol=1)

##Evolutionary rate
plot(NULL, xlim=c(0.5,3.5), ylim=c(0,14), xlab="Dispersal rate",ylab="Evolutionary rate", xaxt='n')
mtext(c("10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:3)
points(c(1:3+gap[1]),apply(ltm_ratio[-6,10:12],2,median)*1000,pch=19, col="orange", cex=1.3)
segments(c(1:3+gap[1]),apply(ltm_ratio[-6,10:12],2,quantile,probs=.025)*1000,
         c(1:3+gap[1]),apply(ltm_ratio[-6,10:12],2,quantile,probs=.975)*1000,lwd=1, col="orange")
segments(c(1:3+gap[1]),apply(ltm_ratio[-6,10:12],2,quantile,probs=.25)*1000,
         c(1:3+gap[1]),apply(ltm_ratio[-6,10:12],2,quantile,probs=.75)*1000,lwd=2, col="orange")

points(c(1:3+gap[2]),apply(ltm_ratio[-6,13:15],2,median)*1000,pch=19, col="red", cex=1.3)
segments(c(1:3+gap[2]),apply(ltm_ratio[-6,13:15],2,quantile,probs=.025)*1000,
         c(1:3+gap[2]),apply(ltm_ratio[-6,13:15],2,quantile,probs=.975)*1000,lwd=1, col="red")
segments(c(1:3+gap[2]),apply(ltm_ratio[-6,13:15],2,quantile,probs=.25)*1000,
         c(1:3+gap[2]),apply(ltm_ratio[-6,13:15],2,quantile,probs=.75)*1000,lwd=2, col="red")

points(c(1:3+gap[3]),apply(ltm_ratio[-6,16:18],2,median)*1000,pch=19, col="dodgerblue", cex=1.3)
segments(c(1:3+gap[3]),apply(ltm_ratio[-6,16:18],2,quantile,probs=.025)*1000,
         c(1:3+gap[3]),apply(ltm_ratio[-6,16:18],2,quantile,probs=.975)*1000,lwd=1, col="dodgerblue")
segments(c(1:3+gap[3]),apply(ltm_ratio[-6,16:18],2,quantile,probs=.25)*1000,
         c(1:3+gap[3]),apply(ltm_ratio[-6,16:18],2,quantile,probs=.75)*1000,lwd=2, col="dodgerblue")

points(c(1:3+gap[4]),apply(ltm_ratio[-6,19:21],2,median)*1000,pch=19, col="darkviolet", cex=1.3)
segments(c(1:3+gap[4]),apply(ltm_ratio[-6,19:21],2,quantile,probs=.025)*1000,
         c(1:3+gap[4]),apply(ltm_ratio[-6,19:21],2,quantile,probs=.975)*1000,lwd=1, col="darkviolet")
segments(c(1:3+gap[4]),apply(ltm_ratio[-6,19:21],2,quantile,probs=.25)*1000,
         c(1:3+gap[4]),apply(ltm_ratio[-6,19:21],2,quantile,probs=.75)*1000,lwd=2, col="darkviolet")

legend("topleft", cex=.8, legend=c("Random div.","Random + Dist","Random + Size","Random +D+S"), col=2:3, border=NA,fill=c("orange","red","dodgerblue","darkviolet"), bty="n", ncol=1)


#################################################################################
### Fig.6 TEMPORAL EVOLUTION OF gG SOME POPULATIONS SCN SPATIAL CONFIGURATION ###
#################################################################################

#graph
par(mfrow=c(3,3))
for (pop in 1:15) {
  geno_1SW_tot<-list()
  for (scn in scn_spat_10){ #define scenarios
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  color_transparent <- function(color) adjustcolor(color, alpha.f = 0.3)
  gap <- c(-.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3,-.2)
  
  plot(NULL, xlim=c(0,50), ylim=c(-.15,.15), xlab="Years",ylab="Growth potential", main=pop, cex.lab=1.5,cex.axis=1.2)
  rect(0,-.15,10,.15, density = NULL, col = color_transparent("grey"), border = NA)
  abline(h=0, lty=2)
  
  for (scn in scn_spat_10) { #define scenarios here too
    median_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    lw<-loess(median_gene~c(1:50),span=0.2)
    lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
    
    #if(scn %in% c(2:4,6:8,10:12) ) {
    if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
      median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
      lines((12:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[10:50,pop,scn]),lty=3)
    }
  }
  #legend("topright", cex=.8, legend=c("No div. - No dispersal","No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Gradual div. - No dispersal","Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[5:8], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - No dispersal","Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[9:12], bty="n", ncol=2)
  legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(1,4,7,10)], bty="n", ncol=1)
  legend("bottomright", cex=.8, legend=c("Philo","Immi"), border=NA,lty=c(1,3), bty="n")
}


#################################################################
### Fig.7 BOXPLOT LOCAL CV SCN SIMPLE / SPATIAL CONFIGURATION ###
#################################################################

##Simple spatial configuration
gap<-c(-.2,0,.2)
plot(NULL, xlim=c(0.5,5.5), ylim=c(0.06,0.12), xlab="Dispersal scenario",ylab="Local CV growth potential", xaxt='n')
mtext(c("Init.", "0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:5)
sd_init_median_nodiv <- apply(cv_init[,c(1:4)],1,mean)
sd_init_median_grad <- apply(cv_init[,c(5:8)],1,mean)
sd_init_median_rand <- apply(cv_init[,c(9:12)],1,mean)
abline(v=1.5,lty=3)
points(c(1:5+gap[1]),c(median(sd_init_median_nodiv),apply(cv_end[,1:4],2,median)),pch=19, cex=1.3)
segments(c(1:5+gap[1]),c(quantile(sd_init_median_nodiv,probs=.025),apply(cv_end[,1:4],2,quantile,probs=.025)),
         c(1:5+gap[1]),c(quantile(sd_init_median_nodiv,probs=.975),apply(cv_end[,1:4],2,quantile,probs=.975)),lwd=1)
segments(c(1:5+gap[1]),c(quantile(sd_init_median_nodiv,probs=.25),apply(cv_end[,1:4],2,quantile,probs=.25)),
         c(1:5+gap[1]),c(quantile(sd_init_median_nodiv,probs=.75),apply(cv_end[,1:4],2,quantile,probs=.75)),lwd=2)

points(c(1:5+gap[2]),c(median(sd_init_median_grad),apply(cv_end[,5:8],2,median)),pch=19, col="forestgreen", cex=1.3)
segments(c(1:5+gap[2]),c(quantile(sd_init_median_grad,probs=.025),apply(cv_end[,5:8],2,quantile,probs=.025)),
         c(1:5+gap[2]),c(quantile(sd_init_median_grad,probs=.975),apply(cv_end[,5:8],2,quantile,probs=.975)),lwd=1, col="forestgreen")
segments(c(1:5+gap[2]),c(quantile(sd_init_median_grad,probs=.25),apply(cv_end[,5:8],2,quantile,probs=.25)),
         c(1:5+gap[2]),c(quantile(sd_init_median_grad,probs=.75),apply(cv_end[,5:8],2,quantile,probs=.75)),lwd=2, col="forestgreen")

points(c(1:5+gap[3]),c(median(sd_init_median_rand),apply(cv_end[,9:12],2,median)),pch=19, col="orange", cex=1.3)
segments(c(1:5+gap[3]),c(quantile(sd_init_median_rand,probs=.025),apply(cv_end[,9:12],2,quantile,probs=.025)),
         c(1:5+gap[3]),c(quantile(sd_init_median_rand,probs=.975),apply(cv_end[,9:12],2,quantile,probs=.975)),lwd=1, col="orange")
segments(c(1:5+gap[3]),c(quantile(sd_init_median_rand,probs=.25),apply(cv_end[,9:12],2,quantile,probs=.25)),
         c(1:5+gap[3]),c(quantile(sd_init_median_rand,probs=.75),apply(cv_end[,9:12],2,quantile,probs=.75)),lwd=2, col="orange")
legend("topright", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=c("black","forestgreen","orange"), bty="n", ncol=1)


##Spatial configuration
gap<-c(-.3,-.1,.1,.3)
plot(NULL, xlim=c(0.5,4.5), ylim=c(0.06,0.12), xlab="Dispersal rate",ylab="Local CV growth potential", xaxt='n')
mtext(c("Init.","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4)
sd_init_median_rand <- apply(cv_init[,c(10:12)],1,mean)
sd_init_median_dist <- apply(cv_init[,c(13:15)],1,mean)
sd_init_median_size <- apply(cv_init[,c(16:18)],1,mean)
sd_init_median_both <- apply(cv_init[,c(19:21)],1,mean)
abline(v=1.5,lty=3)
points(c(1:4+gap[1]),c(median(sd_init_median_rand),apply(cv_end[,10:12],2,median)),pch=19, col="orange", cex=1.3)
segments(c(1:4+gap[1]),c(quantile(sd_init_median_rand,probs=.025),apply(cv_end[,10:12],2,quantile,probs=.025)),
         c(1:4+gap[1]),c(quantile(sd_init_median_rand,probs=.975),apply(cv_end[,10:12],2,quantile,probs=.975)),lwd=1, col="orange")
segments(c(1:4+gap[1]),c(quantile(sd_init_median_rand,probs=.25),apply(cv_end[,10:12],2,quantile,probs=.25)),
         c(1:4+gap[1]),c(quantile(sd_init_median_rand,probs=.75),apply(cv_end[,10:12],2,quantile,probs=.75)),lwd=2, col="orange")

points(c(1:4+gap[2]),c(median(sd_init_median_dist),apply(cv_end[,13:15],2,median)),pch=19, col="red", cex=1.3)
segments(c(1:4+gap[2]),c(quantile(sd_init_median_dist,probs=.025),apply(cv_end[,13:15],2,quantile,probs=.025)),
         c(1:4+gap[2]),c(quantile(sd_init_median_dist,probs=.975),apply(cv_end[,13:15],2,quantile,probs=.975)),lwd=1, col="red")
segments(c(1:4+gap[2]),c(quantile(sd_init_median_dist,probs=.25),apply(cv_end[,13:15],2,quantile,probs=.25)),
         c(1:4+gap[2]),c(quantile(sd_init_median_dist,probs=.75),apply(cv_end[,13:15],2,quantile,probs=.75)),lwd=2, col="red")

points(c(1:4+gap[3]),c(median(sd_init_median_size),apply(cv_end[,16:18],2,median)),pch=19, col="dodgerblue", cex=1.3)
segments(c(1:4+gap[3]),c(quantile(sd_init_median_size,probs=.025),apply(cv_end[,16:18],2,quantile,probs=.025)),
         c(1:4+gap[3]),c(quantile(sd_init_median_size,probs=.975),apply(cv_end[,16:18],2,quantile,probs=.975)),lwd=1, col="dodgerblue")
segments(c(1:4+gap[3]),c(quantile(sd_init_median_size,probs=.25),apply(cv_end[,16:18],2,quantile,probs=.25)),
         c(1:4+gap[3]),c(quantile(sd_init_median_size,probs=.75),apply(cv_end[,16:18],2,quantile,probs=.75)),lwd=2, col="dodgerblue")

points(c(1:4+gap[4]),c(median(sd_init_median_both),apply(cv_end[,19:21],2,median)),pch=19, col="darkviolet", cex=1.3)
segments(c(1:4+gap[4]),c(quantile(sd_init_median_both,probs=.025),apply(cv_end[,19:21],2,quantile,probs=.025)),
         c(1:4+gap[4]),c(quantile(sd_init_median_both,probs=.975),apply(cv_end[,19:21],2,quantile,probs=.975)),lwd=1, col="darkviolet")
segments(c(1:4+gap[4]),c(quantile(sd_init_median_both,probs=.25),apply(cv_end[,19:21],2,quantile,probs=.25)),
         c(1:4+gap[4]),c(quantile(sd_init_median_both,probs=.75),apply(cv_end[,19:21],2,quantile,probs=.75)),lwd=2, col="darkviolet")

legend("topright", cex=.8, legend=c("Random div.","Random + Dist","Random + Size","Random +D+S"), col=2:3, border=NA,fill=c("orange","red","dodgerblue","darkviolet"), bty="n", ncol=1)




################################## SUPPLEMENTARY MATERIALS ####################################

##################################################################
### Fig.S1 TEMPORAL EVOLUTION gG ALL POPULATIONS SCN DIVERSITY ###
##################################################################

#graph
par(mfrow=c(3,3))
for (pop in 1:15) {
  geno_1SW_tot<-list()
  for (scn in scn_all_div){ #define scenarios
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  color_transparent <- function(color) adjustcolor(color, alpha.f = 0.3)
  gap <- c(-.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3,-.2)
  
  plot(NULL, xlim=c(0,50), ylim=c(-.15,.15), xlab="Years",ylab="Growth potential", main=pop, cex.lab=1.5,cex.axis=1.2)
  rect(0,-.15,10,.15, density = NULL, col = color_transparent("grey"), border = NA)
  abline(h=0, lty=2)
  
  for (scn in scn_all_div) { #define scenarios here too
    median_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    lw<-loess(median_gene~c(1:50),span=0.2)
    lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
    
    if(scn %in% c(2:4,6:8,10:12) ) {
    #if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
      median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
      lines((12:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[10:50,pop,scn]),lty=3)
    }
  }
  #legend("topright", cex=.8, legend=c("No div. - No dispersal","No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Gradual div. - No dispersal","Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[5:8], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - No dispersal","Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[9:12], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(1,4,7,10)], bty="n", ncol=1)
  legend("bottomright", cex=.8, legend=c("Philo","Immi"), border=NA,lty=c(1,3), bty="n")
}

##############################################################################
### Fig.S2 TEMPORAL EVOLUTION gG ALL POPULATIONS SCN SPATIAL CONFIGURATION ###
##############################################################################

#graph
par(mfrow=c(3,3))
for (pop in 1:15) {
  geno_1SW_tot<-list()
  for (scn in scn_spat_10){ #define scenarios
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  color_transparent <- function(color) adjustcolor(color, alpha.f = 0.3)
  gap <- c(-.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3,-.2)
  
  plot(NULL, xlim=c(0,50), ylim=c(-.15,.15), xlab="Years",ylab="Growth potential", main=pop, cex.lab=1.5,cex.axis=1.2)
  rect(0,-.15,10,.15, density = NULL, col = color_transparent("grey"), border = NA)
  abline(h=0, lty=2)
  
  for (scn in scn_spat_10) { #define scenarios here too
    median_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    lw<-loess(median_gene~c(1:50),span=0.2)
    lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
    
    #if(scn %in% c(2:4,6:8,10:12) ) {
      if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
      median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
      lines((12:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[10:50,pop,scn]),lty=3)
    }
  }
  #legend("topright", cex=.8, legend=c("No div. - No dispersal","No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Gradual div. - No dispersal","Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[5:8], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - No dispersal","Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[9:12], bty="n", ncol=2)
  legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(1,4,7,10)], bty="n", ncol=1)
  legend("bottomright", cex=.8, legend=c("Philo","Immi"), border=NA,lty=c(1,3), bty="n")
}


##################################################################
### Fig.S3 TEMPORAL EVOLUTION pG ALL POPULATIONS SCN DIVERSITY ###
##################################################################

#par(mfrow=c(2,3))
for (pop in 1:15) {
  geno_1SW_tot<-list()
  for (scn in scn_all_div){ #define scenarios
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  color_transparent <- function(color) adjustcolor(color, alpha.f = 0.3)
  gap <- c(-.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3)
  
  
  plot(NULL, xlim=c(0,50), ylim=c(-0.15,0.15), xlab="Years",ylab="Growth potential phenotype", main=pop, cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(1,1.5), xlab="Years",ylab="Adult gFmid1", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(33,47), xlab="Years",ylab="Adult gFmid2", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(30,50), xlab="Years",ylab="Adult gFmid3", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(70,100), xlab="Years",ylab="Adult gFmid4", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(0.2,.8), xlab="Years",ylab="Adult gNeutral", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  rect(0,-.2,10,.2, density = NULL, col = color_transparent("grey"), border = NA)
  abline(h=0, lty=2)
  
  for (scn in scn_all_div) { #define scenarios here too
    median_gene<-aggregate(geno_1SW_tot[[scn]]$pG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid1[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid2[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid3[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid4[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gNeutral[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    lw<-loess(median_gene~c(1:50),span=0.2)
    lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=2)
    
    if(scn %in% c(2:4,6:8,10:12) ) {
    #if(scn >1) { #immigrants
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$pG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid1[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid2[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid3[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid4[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gNeutral[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
      lines((12:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[10:50,pop,scn]),lty=3)
    }
  }
  #legend("topright", cex=.8, legend=c("No div. - No dispersal","No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Gradual div. - No dispersal","Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - No dispersal","Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[1:4], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(1,4,7,10)], bty="n", ncol=1)
  legend("bottomright", cex=.8, legend=c("Philo","Immi"), border=NA,lty=c(1,3), bty="n")
}

##############################################################################
### Fig.S4 TEMPORAL EVOLUTION pG ALL POPULATIONS SCN SPATIAL CONFIGURATION ###
##############################################################################

#par(mfrow=c(2,3))
for (pop in 1:15) {
  geno_1SW_tot<-list()
  for (scn in scn_spat_10){ #define scenarios
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  color_transparent <- function(color) adjustcolor(color, alpha.f = 0.3)
  gap <- c(-.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3, -.2,-.1,0, .1, .2, .3)
  
  
  plot(NULL, xlim=c(0,50), ylim=c(-0.15,0.15), xlab="Years",ylab="Growth potential phenotype", main=pop, cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(1,1.5), xlab="Years",ylab="Adult gFmid1", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(33,47), xlab="Years",ylab="Adult gFmid2", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(30,50), xlab="Years",ylab="Adult gFmid3", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(70,100), xlab="Years",ylab="Adult gFmid4", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  #plot(NULL, xlim=c(0,50), ylim=c(0.2,.8), xlab="Years",ylab="Adult gNeutral", main=paste0(pops[pop] ), cex.lab=1,cex.axis=1)
  rect(0,-.2,10,.2, density = NULL, col = color_transparent("grey"), border = NA)
  abline(h=0, lty=2)
  
  for (scn in scn_spat_10) { #define scenarios here too
    median_gene<-aggregate(geno_1SW_tot[[scn]]$pG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid1[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid2[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid3[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gFmid4[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    #median_gene<-aggregate(geno_1SW_tot[[scn]]$gNeutral[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:50,2]
    lw<-loess(median_gene~c(1:50),span=0.2)
    lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=2)
    
    #if(scn %in% c(2:4,6:8,10:12) ) {
    if(scn >9) { #immigrants
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$pG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid1[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid2[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid3[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      #median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gFmid4[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gNeutral[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:50,2]
      lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
      lines((12:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[10:50,pop,scn]),lty=3)
    }
  }
  #legend("topright", cex=.8, legend=c("No div. - No dispersal","No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Gradual div. - No dispersal","Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[1:5], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div. - No dispersal","Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[1:4], bty="n", ncol=2)
  legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(1,4,7,10)], bty="n", ncol=1)
  legend("bottomright", cex=.8, legend=c("Philo","Immi"), border=NA,lty=c(1,3), bty="n")
}


####################################################################
### Fig.S5 RELATIONSHIPS BETWEEN METRICS AND IMMIGRANTS FEATURES ###
####################################################################

## METAPOPULATION TRAIT MISMATCH ##

trait_mismatch4<-diff_trait2<-array(,dim=c(nYears+nInit, npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (pop in 1:npop) {
    geno_1SW<-NULL
    for (simul in 1:nSIMUL) { #keep immigrants from all simulations to have enough of them
      geno_1SW <- rbind(geno_1SW,res_1SW_scn[[scn]][[simul]][[pop]][res_1SW_scn[[scn]][[simul]][[pop]]$CollecID!=pop,])
    }
    for (t in 1:(nYears+nInit)) {
      geno_bis<-geno_1SW$gG[geno_1SW$years==t]
      
      diff_trait <- (optimum[pop]-mean(geno_bis)) #abs
      
      #trait_mismatch4[t,pop,scn]<- abs(diff_trait*P.im.median[t,pop,scn])
      diff_trait2[t,pop,scn]<- abs(diff_trait)
    }
  }
}

#Local trait mismatch function metapopulation trait mismatch
##with ellipses
library(reshape2)
ltm_immi <- array(,dim=c(npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  ltm_immi[,scn] <- apply(diff_trait2[45:50,,scn],2,mean)
}
data_test <- cbind(melt(ltm_end),melt(ltm_immi))
colnames(data_test) <- c("pop","scn","ltm","pop2","scn2","ltm_immi")
data_test$scn_type <- NA
for (i in 1:nrow(data_test)) {
  if (data_test$scn[i] %in% c(2:4)) {
    data_test$scn_type[i] <- "nodiv"
  }
  if (data_test$scn[i] %in% c(6:8)) {
    data_test$scn_type[i] <- "gradual"
  }
  if (data_test$scn[i] %in% c(10:12)) {
    data_test$scn_type[i] <- "random"
  }
  if (data_test$scn[i] %in% c(13:15)) {
    data_test$scn_type[i] <- "random_dist"
  }
  if (data_test$scn[i] %in% c(16:18)) {
    data_test$scn_type[i] <- "random_size"
  }
  if (data_test$scn[i] %in% c(19:21)) {
    data_test$scn_type[i] <- "random_d_s"
  }
}

library(ggplot2)
ggplot(data_test, aes(x = ltm_immi, y = ltm, color = scn_type)) +
  geom_point() +
  stat_ellipse(geom = "polygon",
               aes(fill = scn_type), 
               alpha = 0.1) +
  theme_bw() +
  scale_fill_manual(values=c("black","orange","darkviolet","red","dodgerblue","forestgreen"))+
  scale_color_manual(values=c("black","orange","darkviolet","red","dodgerblue","forestgreen"))+
  xlab("Metapopulation trait mismatch") +
  ylab("Local trait mismatch") +
  theme(axis.text = element_text(size=11), axis.title= element_text(size=13))


#Evolutionary rate function of proportion of immigrants
##with ellipses
library(reshape2)
prop_immi <- array(,dim=c(npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  prop_immi[,scn] <- apply(P.im.median[10:50,,scn],2,mean)
}
data_test <- cbind(melt(ltm_ratio),melt(prop_immi))
colnames(data_test) <- c("pop","scn","ltm_ratio","pop2","scn2","prop_immi")
data_test$ltm_ratio<-data_test$ltm_ratio*1000
data_test<-subset(data_test, data_test$pop !=6)
data_test<-subset(data_test, data_test$pop !=7)
data_test<-subset(data_test, data_test$pop !=8)
data_test$scn_type <- NA
for (i in 1:nrow(data_test)) {
  if (data_test$scn[i] %in% c(2:4)) {
    data_test$scn_type[i] <- "nodiv"
  }
  if (data_test$scn[i] %in% c(6:8)) {
    data_test$scn_type[i] <- "gradual"
  }
  if (data_test$scn[i] %in% c(10:12)) {
    data_test$scn_type[i] <- "random"
  }
  if (data_test$scn[i] %in% c(13:15)) {
    data_test$scn_type[i] <- "random_dist"
  }
  if (data_test$scn[i] %in% c(16:18)) {
    data_test$scn_type[i] <- "random_size"
  }
  if (data_test$scn[i] %in% c(19:21)) {
    data_test$scn_type[i] <- "random_d_s"
  }
}

try<- subset(data_test, !is.na(data_test$scn_type))
library(ggplot2)
ggplot(try, aes(x = prop_immi, y = ltm_ratio, color = scn_type)) +
  geom_point() +
  stat_ellipse(geom = "polygon",
               aes(fill = scn_type), 
               alpha = 0.1) +
  theme_bw() +
  scale_fill_manual(values=c("black","orange","darkviolet","red","dodgerblue","forestgreen"))+
  scale_color_manual(values=c("black","orange","darkviolet","red","dodgerblue","forestgreen"))+
  xlab("Proportion of immigrants") +
  ylab("Evolutionary rate") +
  theme(axis.text = element_text(size=11), axis.title= element_text(size=13))


##################################################################
### Fig.S6 TEMPORAL EVOLUTION OF POPULATION SIZE SCN DIVERSITY ###
##################################################################

par(mfrow=c(3,3))
for (pop in 1:15){
  Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(EXPE)))
  for (scn in 1:length(EXPE)){
    for (simul in 1:nSIMUL){
      Npops[,simul,scn] <- nReturns_scn[[scn]][[1]][[simul]][,pop] #nb returns metapop per year, simu and scenario
    }
  }
  
  plot(NULL, xlim=c(0,(nYears+nInit)), ylim=c(0,200), xlab="Years",ylab="Population abundance",main=pop, cex.lab=1.5,cex.axis=1.2) #500
  rect(0,0,10,200, density = NULL, col = color_transparent("grey"), border = NA)
  
  for (scn in scn_all_div) { #choice of scenarios to plot
      median_npops<- apply(Npops[,,scn],1,median)
      lw<-loess(median_npops~c(1:50),span=0.2)
      lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=2)
  }
  #legend("topright", cex=.8, legend=c("Gradual div.","Random div.","Random + Dist","Size"), col=2:3, border=NA,fill=color_scn[1:4], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div.-0%","Random div.-10%","Random div.-20%","Random div.-30%"), col=2:3, border=NA,fill=color_scn[scn_rand], bty="n", ncol=2)
}

##############################################################################
### Fig.S7 TEMPORAL EVOLUTION OF POPULATION SIZE SCN SPATIAL CONFIGURATION ###
##############################################################################

par(mfrow=c(3,3))
for (pop in 1:15){
  Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(EXPE)))
  for (scn in 1:length(EXPE)){
    for (simul in 1:nSIMUL){
      Npops[,simul,scn] <- nReturns_scn[[scn]][[1]][[simul]][,pop] #nb returns metapop per year, simu and scenario
    }
  }
  
  plot(NULL, xlim=c(0,(nYears+nInit)), ylim=c(0,200), xlab="Years",ylab="Population abundance",main=pop, cex.lab=1.5,cex.axis=1.2) #500
  rect(0,0,10,200, density = NULL, col = color_transparent("grey"), border = NA)
  
  for (scn in scn_spat_10) { #choice of scenarios to plot
    median_npops<- apply(Npops[,,scn],1,median)
    lw<-loess(median_npops~c(1:50),span=0.2)
    lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=2)
  }
  #legend("topright", cex=.8, legend=c("Gradual div.","Random div.","Random + Dist","Size"), col=2:3, border=NA,fill=color_scn[1:4], bty="n", ncol=2)
  #legend("topright", cex=.8, legend=c("Random div.-0%","Random div.-10%","Random div.-20%","Random div.-30%"), col=2:3, border=NA,fill=color_scn[scn_rand], bty="n", ncol=2)
}


###############################################
### Fig.S8 METAPOPULATION SIZE / CV ALL SCN ###
###############################################

nMetapop_size_5 <-nMetapop_cv<-nMetapop_var <-nMetapop_mean<- array(,dim=c(nSIMUL,length(EXPE)))
Npops_mean<-Npops_var<-array(,dim=c(npop,nSIMUL,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    nMetapop_size_5[simul,scn]<-median(rowSums(nReturns_scn[[scn]][[1]][[simul]])[45:50])
    nMetapop_cv[simul,scn]<-cv(rowSums(nReturns_scn[[scn]][[1]][[simul]])[10:50])
    nMetapop_var[simul,scn]<-var(rowSums(nReturns_scn[[scn]][[1]][[simul]])[10:50])
    nMetapop_mean[simul,scn]<-mean(rowSums(nReturns_scn[[scn]][[1]][[simul]])[10:50])
  }
}

#Scn diversity simple spatial configuration
gap <- c(-.15,0,.15)
plot(NULL, xlim=c(0.5,4.5), ylim=c(0,3000), xlab="Dispersal rate",ylab="Metapopulation size", xaxt='n')
mtext(c("0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4)
for (scn in 1:4) {
  points(c(scn+gap[1]),median(nMetapop_size_5[,scn]),pch=19, cex=1.3)
  segments(c(scn+gap[1]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn+gap[1]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1)
  segments(c(scn+gap[1]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn+gap[1]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2)
}
for (scn in 5:8) {
  points(c(scn-4+gap[2]),median(nMetapop_size_5[,scn]),pch=19, col="forestgreen", cex=1.3)
  segments(c(scn-4+gap[2]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn-4+gap[2]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1, col="forestgreen")
  segments(c(scn-4+gap[2]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn-4+gap[2]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2, col="forestgreen")
}
for (scn in 9:12) {
  points(c(scn-8+gap[3]),median(nMetapop_size_5[,scn]),pch=19, col="orange", cex=1.3)
  segments(c(scn-8+gap[3]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn-8+gap[3]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1, col="orange")
  segments(c(scn-8+gap[3]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn-8+gap[3]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2, col="orange")
}
legend("topleft", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=c("black","forestgreen","orange"), bty="n", ncol=1)


plot(NULL, xlim=c(0.5,4.5), ylim=c(0,0.15), xlab="Dispersal rate",ylab="Metapopulation CV", xaxt='n')
mtext(c("0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4)
for (scn in 1:4) {
  points(c(scn+gap[1]),median(nMetapop_cv[,scn]),pch=19, cex=1.3)
  segments(c(scn+gap[1]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn+gap[1]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1)
  segments(c(scn+gap[1]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn+gap[1]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2)
}
for (scn in 5:8) {
  points(c(scn-4+gap[2]),median(nMetapop_cv[,scn]),pch=19, col="forestgreen", cex=1.3)
  segments(c(scn-4+gap[2]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn-4+gap[2]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1, col="forestgreen")
  segments(c(scn-4+gap[2]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn-4+gap[2]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2, col="forestgreen")
}
for (scn in 9:12) {
  points(c(scn-8+gap[3]),median(nMetapop_cv[,scn]),pch=19, col="orange", cex=1.3)
  segments(c(scn-8+gap[3]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn-8+gap[3]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1, col="orange")
  segments(c(scn-8+gap[3]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn-8+gap[3]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2, col="orange")
}
legend("bottomleft", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=c("black","forestgreen","orange"), bty="n", ncol=1)


#scn spatial configuration
gap<-c(-0.15,-.05,0.05,.15)
plot(NULL, xlim=c(0.5,3.5), ylim=c(0,3000), xlab="Dispersal rate",ylab="Metapopulation size", xaxt='n')
mtext(c("10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:3)
for (scn in c(10,11,12)) { #rand
  points(c(scn-9+gap[1]),median(nMetapop_size_5[,scn]),pch=19, col="orange", cex=1.3)
  segments(c(scn-9+gap[1]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn-9+gap[1]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1, col="orange")
  segments(c(scn-9+gap[1]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn-9+gap[1]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2, col="orange")
}
for (scn in c(13,14,15)) { #dist
  points(c(scn-12+gap[2]),median(nMetapop_size_5[,scn]),pch=19, col="red", cex=1.3)
  segments(c(scn-12+gap[2]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn-12+gap[2]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1, col="red")
  segments(c(scn-12+gap[2]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn-12+gap[2]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2, col="red")
}
for (scn in c(16,17,18)) { #size
  points(c(scn-15+gap[3]),median(nMetapop_size_5[,scn]),pch=19, col="dodgerblue", cex=1.3)
  segments(c(scn-15+gap[3]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn-15+gap[3]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1, col="dodgerblue")
  segments(c(scn-15+gap[3]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn-15+gap[3]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2, col="dodgerblue")
}
for (scn in c(19,20,21)) { #both
  points(c(scn-18+gap[4]),median(nMetapop_size_5[,scn]),pch=19, col="darkviolet", cex=1.3)
  segments(c(scn-18+gap[4]),quantile(nMetapop_size_5[,scn],probs=.025),
           c(scn-18+gap[4]),quantile(nMetapop_size_5[,scn],probs=.975),lwd=1, col="darkviolet")
  segments(c(scn-18+gap[4]),quantile(nMetapop_size_5[,scn],probs=.25),
           c(scn-18+gap[4]),quantile(nMetapop_size_5[,scn],probs=.75),lwd=2, col="darkviolet")
}
legend("bottomright", cex=.8, legend=c("Random div.","Random + Dist","Random + Size","Random +D+S"), col=2:3, border=NA,fill=c("orange","coral1","darkcyan","hotpink4"), bty="n", ncol=1)


plot(NULL, xlim=c(0.5,3.5), ylim=c(0,0.15), xlab="Dispersal rate",ylab="Metapopulation CV", xaxt='n')
mtext(c("10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:3)
for (scn in c(10,11,12)) { #rand
  points(c(scn-9+gap[1]),median(nMetapop_cv[,scn]),pch=19, col="orange", cex=1.3)
  segments(c(scn-9+gap[1]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn-9+gap[1]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1, col="orange")
  segments(c(scn-9+gap[1]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn-9+gap[1]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2, col="orange")
}
for (scn in c(13,14,15)) { #dist
  points(c(scn-12+gap[2]),median(nMetapop_cv[,scn]),pch=19, col="red", cex=1.3)
  segments(c(scn-12+gap[2]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn-12+gap[2]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1, col="red")
  segments(c(scn-12+gap[2]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn-12+gap[2]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2, col="red")
}
for (scn in c(16,17,18)) { #size
  points(c(scn-15+gap[3]),median(nMetapop_cv[,scn]),pch=19, col="dodgerblue", cex=1.3)
  segments(c(scn-15+gap[3]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn-15+gap[3]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1, col="dodgerblue")
  segments(c(scn-15+gap[3]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn-15+gap[3]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2, col="dodgerblue")
}
for (scn in c(19,20,21)) { #both
  points(c(scn-18+gap[4]),median(nMetapop_cv[,scn]),pch=19, col="darkviolet", cex=1.3)
  segments(c(scn-18+gap[4]),quantile(nMetapop_cv[,scn],probs=.025),
           c(scn-18+gap[4]),quantile(nMetapop_cv[,scn],probs=.975),lwd=1, col="darkviolet")
  segments(c(scn-18+gap[4]),quantile(nMetapop_cv[,scn],probs=.25),
           c(scn-18+gap[4]),quantile(nMetapop_cv[,scn],probs=.75),lwd=2, col="darkviolet")
}
legend("bottomright", cex=.8, legend=c("Random div.","Random + Dist","Random + Size","Random +D+S"), col=2:3, border=NA,fill=c("orange","coral1","darkcyan","hotpink4"), bty="n", ncol=1)


###############################################################
### Fig.S9 RATIO IMMIGRANTS/EMIGRANTS BY POPULATION ALL SCN ###
###############################################################

I <- P.im <- array(,dim=c(nSIMUL,npop, length(EXPE)))
ratio2 <- array(,dim=c((nYears+nInit),npop,nSIMUL, length(EXPE)))
for (scn in c(2:4,6:8,10:21)){
  for (simul in 1:nSIMUL){
    ratio<-Mig_scn[[scn]][[1]][[simul]]$NIm/Mig_scn[[scn]][[1]][[simul]]$NEm
    ratio[which(ratio=="Inf")]<-NA
    I[simul, ,scn] <- apply(ratio[45:50,],2,median, na.rm=T) #5last years #mean , na.rm=T
  }
}

gap<-c(-0.3,-0.15,  0.0,  0.15, -0.3, 0,  0.0,  0.15, -0.3, 0.15,  0.0,  0.15, 0,  0.0,  0.15, 0.15,  0.0,  0.15, 0.3,  0.0,  0.15)
plot(NULL, xlim=c(0,16), ylim=c(0,9), xlab="Population",ylab="Ratio Immigrants/Emigrants", xaxt='n')
mtext(1:15, side = 1, line = 1, outer = FALSE, at = 1:15)
abline(h=1,col="black",lty=2)
color_scn_here=c("","black","","","","forestgreen","","","","orange")
for (scn in c(2,6,10)){ #scn_spat_10
  for (pop in 1:npop) {
    points(pop+gap[scn], median(I[,pop,scn],na.rm=TRUE), pch=20, cex=1.8,col=color_scn_here[scn]);
    segments(pop+gap[scn], quantile(I[,pop,scn],probs=.025,na.rm=TRUE),pop+gap[scn],quantile(I[,pop,scn],probs=.975,na.rm=TRUE),col=color_scn_here[scn])
    segments(pop+gap[scn], quantile(I[,pop,scn],probs=.25,na.rm=TRUE),pop+gap[scn],quantile(I[,pop,scn],probs=.75,na.rm=TRUE),col=color_scn_here[scn], lwd=2)
  }
}

#legend("topright", cex=.8, legend=c("No div. - 10%","No div. - 20%", "No div. - 30%"), col=2:3, border=NA,fill=color_scn[2:5], bty="n", ncol=2)
#legend("topright", cex=.8, legend=c("Gradual div. - 10%","Gradual div. - 20%", "Gradual div. - 30%"), col=2:3, border=NA,fill=color_scn[6:8], bty="n", ncol=2)
#legend("topright", cex=.8, legend=c("Random div. - 10% dispersal", "Random div. - 20% dispersal", "Random div. - 30% dispersal"), col=2:3, border=NA,fill=color_scn[10:12], bty="n", ncol=2)
#legend("topright", cex=.8, legend=c("Random div. - 10%","Random + Dist - 10%","Random + Size - 10%", "Random + D+S - 10%"), col=2:3, border=NA,fill=color_scn[c(10,13,16,19)], bty="n", ncol=1)


################################
### Fig.S10 NETWORKS ALL SCN ###
################################

###########  Migrants flows

a <- array(,dim=c(npop,npop,nSIMUL, length(EXPE)))
for (scn in c(2:4,6:8,10:21)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      a[pop,,simul,scn]<-apply(Mig_scn[[scn]][[1]][[simul]]$Im[[pop]][46:50,],2,mean)
      for (pop2 in 1:npop) {
        if (a[pop,pop2,simul,scn]=="NaN" || is.na(a[pop,pop2,simul,scn])) {
          a[pop,pop2,simul,scn]=0
        }
      }
    }
  }
}
arr<-mflow<-list()
for (scn in 1:length(EXPE)) {
  arr[[scn]] <- array( unlist(a[,,,scn]) , c(15,15,100) )
  mflow[[scn]]<-apply( arr[[scn]] , 1:2 , mean )#mean simulations
  colnames(mflow[[scn]])<-pops
  rownames(mflow[[scn]])<-pops
}


mflow2<-melt(mflow[[2]]) #for absolute values #SCENARIO CHOICE #2,6,10,13,16,19
mflow3<-mflow2[c(2,1,3)] #for absolute values

mflow3[,1]<-as.factor(mflow3[,1])
mflow3[,2]<-as.factor(mflow3[,2])

########### Populations size
Parpop<-array(dim=c(nSIMUL, npop, length(EXPE)))
Parpop2<-array(dim=c(npop, length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      Parpop[simul, ,scn] <- colMeans(nReturns_scn[[scn]][[1]][[simul]][46:50,], na.rm=T) #5last years
    }
  }
  Parpop2[,scn]<-colMeans(Parpop[,,scn])
}
Type[7]<-"source"
Type[c(1,2,4,7,8,9)]<-c("source","sink","sink","source","source","sink")

nodes<-data.frame(id=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15"),
                  pop=pops,
                  type=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"), #control
                  size=Parpop2[,2]#choice of scenarios #2,6,10,13,16,19
)

############### Network
mflow3$Var2=as.character(mflow3$Var2)
mflow3$Var1=as.character(mflow3$Var1)

mflow3$Var2[which(mflow3$Var2=="Leff")]<-"p1" ;mflow3$Var1[which(mflow3$Var1=="Leff")]<-"p1"
mflow3$Var2[which(mflow3$Var2=="Trieux")]<-"p2" ;mflow3$Var1[which(mflow3$Var1=="Trieux")]<-"p2"
mflow3$Var2[which(mflow3$Var2=="Jaudy")]<-"p3" ;mflow3$Var1[which(mflow3$Var1=="Jaudy")]<-"p3"
mflow3$Var2[which(mflow3$Var2=="Leguer")]<-"p4" ;mflow3$Var1[which(mflow3$Var1=="Leguer")]<-"p4"
mflow3$Var2[which(mflow3$Var2=="Yar")]<-"p5" ;mflow3$Var1[which(mflow3$Var1=="Yar")]<-"p5"
mflow3$Var2[which(mflow3$Var2=="Douron")]<-"p6" ;mflow3$Var1[which(mflow3$Var1=="Douron")]<-"p6"
mflow3$Var2[which(mflow3$Var2=="Penze")]<-"p7" ;mflow3$Var1[which(mflow3$Var1=="Penze")]<-"p7"
mflow3$Var2[which(mflow3$Var2=="Elorn")]<-"p8" ;mflow3$Var1[which(mflow3$Var1=="Elorn")]<-"p8"
mflow3$Var2[which(mflow3$Var2=="Aulne")]<-"p9" ;mflow3$Var1[which(mflow3$Var1=="Aulne")]<-"p9"
mflow3$Var2[which(mflow3$Var2=="Goyen")]<-"p10" ;mflow3$Var1[which(mflow3$Var1=="Goyen")]<-"p10"
mflow3$Var2[which(mflow3$Var2=="Odet")]<-"p11" ;mflow3$Var1[which(mflow3$Var1=="Odet")]<-"p11"
mflow3$Var2[which(mflow3$Var2=="Aven")]<-"p12" ;mflow3$Var1[which(mflow3$Var1=="Aven")]<-"p12"
mflow3$Var2[which(mflow3$Var2=="Laita")]<-"p13" ;mflow3$Var1[which(mflow3$Var1=="Laita")]<-"p13"
mflow3$Var2[which(mflow3$Var2=="Scorff")]<-"p14" ;mflow3$Var1[which(mflow3$Var1=="Scorff")]<-"p14"
mflow3$Var2[which(mflow3$Var2=="Blavet")]<-"p15" ;mflow3$Var1[which(mflow3$Var1=="Blavet")]<-"p15"

links<-mflow3
colnames(links)<-c("from","to","weigth")
links<-links[-which(links$weigth==0),]
vis.nodes<-nodes
vis.links<-links

vis.nodes$shape  <- "dot"  
vis.nodes$shadow <- TRUE # Nodes will drop shadow
#vis.nodes$title  <- nodes$pop # Text on click
vis.nodes$label  <- nodes$pop # Node label
vis.nodes$size   <- nodes$size/5 # Node size
vis.nodes$borderWidth <- 2 # Node border width

#for scn random diversity
gradient_color_random <- c("#CCE1E2", "#E2BAD5","#D17EB8","#8BC9CB", #new colors
                           "#ADD5D7","#E9EEEE","#EDE2E9","#CA69AF",
                           "#D793C1","#009B9F","#DDA7CB","#00A6AA",
                           "#E8CEDF","#62BDC0","#1EB2B5")
#for scn gradual diversity
gradient_color_gradual <- c("#009B9F", "#00A6AA", "#1EB2B5", "#62BDC0",
                            "#8BC9CB", "#ADD5D7", "#CCE1E2", "#E9EEEE",
                            "#EDE2E9", "#E8CEDF", "#E2BAD5", "#DDA7CB",
                            "#D793C1", "#D17EB8", "#CA69AF")
#for scn no diversity
gradient_color_nodiv <- rep(gradient_color_gradual[1],15)

#here put gradient color of the scn
vis.nodes$color.background <- gradient_color_nodiv #rev(gradient_color_func(15))[as.numeric(as.character(nodes$type))]
vis.nodes$color.border <- "black" #rev(gradient_color_func(15))[as.numeric(as.character(nodes$type))]#"orange"
vis.nodes$borderWidth <- .7 #rev(gradient_color_func(15))[as.numeric(as.character(nodes$type))]#"orange"
vis.nodes$color.highlight.background <- "red"
vis.nodes$color.highlight.border <- "gold"
tryagain<-links$weigth*2 #+1 for absolute values #*10 for proportion
vis.links$width <- tryagain # line width
vis.links$arrows <- "to" # arrows: 'from', 'to', or 'middle'
vis.links$smooth <- TRUE    # should the edges be curved?
vis.links$shadow <- FALSE    # edge shadow

#nodes (pop) geographical coordinates
#if complex distances
lat<-c(48.701791,48.677231, 48.714110, 48.651108, 48.646443, 48.636448, 48.600963, 48.479118, 48.204225, 48.040924, 48.002719, 47.869874, 47.870947, 47.859305, 47.834972)
lon<-c(-3.056085, -3.157408 , -3.262510, -3.420931, -3.577629, -3.659678, -3.937004, -4.191071, -4.050434,-4.478585, -4.111829, -3.726740,  -3.545196, -3.400563, -3.207688)
#plot(lon, lat)

#if equal distances
lat2<-seq(max(lat),min(lat),by=-(max(lat)-min(lat))/14)
lon2<-c(seq(max(lon),min(lon),by=-(max(lon)-min(lon))/7),seq(min(lon),max(lon),by=(max(lon)-min(lon))/7))
lon2<-lon2[-9]
#plot(lon2, lat2)

#here put coordinates wanted
vis.nodes$x<- lon2*1000
vis.nodes$y <- -lat2*1000

a<-visNetwork(vis.nodes, vis.links)
a<-visEdges(a, arrows=list(to=list(enable=T, scaleFactor=1.5)),color = list(color = "black", highlight = "red"), smooth = list(enabled = TRUE, type = "diagonalCross"))
a<-visNodes(a,fixed = TRUE,physics=T, font=list(color="black", size=0))
a


####################################
### Fig.S11 REPRODUCTIVE SUCCESS ###
####################################

##focus on scn dispersal 10% random diversity simple spatial configuration

load(paste0("results/FITNESS",30331211,"_50.RData"))

female_simul<-male_simul<-bothsex<-NULL
selection_intensity_fem<-array(,dim=c(45,npop,nSIMUL))
selection_intensity_mal<-array(,dim=c(45,npop,nSIMUL))
selection_intensity_bothsex<-array(,dim=c(45,npop,nSIMUL))

for (simul in 1:50) {
  cat("Simulation :",simul," / ")
  test<-NULL
  for (pop in 1:15) {
    adults<- res_adult[[simul]][[pop]][,c("ID","years","Female","CollecID","pG","motherID","fatherID")]
    adults<-adults[-which(duplicated(adults$ID)),]
    test <- rbind(test, cbind(adults,pop)) #get adults of pops
  }
  mothers<-melt(table(test$motherID)) #nb of offspring per mother
  colnames(mothers)<-c("ID","Count")
  female<-merge(test[which(test$Female==1),], mothers, by="ID", all.x=T) #associate this mother with its information (its current pop, its origin)
  
  fathers<-melt(table(test$fatherID)) #nb of offspring per father
  colnames(fathers)<-c("ID","Count")
  male<-merge(test[which(test$Female==0),], fathers, by="ID", all.x=T)
  
  female$fitness <- na_replace(female$Count,0)
  male$fitness <- na_replace(male$Count,0)
  
  female_simul<-rbind(female_simul, female) #pool simul
  male_simul<-rbind(male_simul, male)
  
  bothsex <- rbind(bothsex,female, male)
}

bothsex_10 <- bothsex

female_10 <- female_simul
male_10 <- male_simul



## graph along time Mean LRS
#comparison random dispersal philo / immigrants

par(mfrow = c(3,3))#it goes c(bottom, left, top, right)
for (popo in 1:15) {
  
  #bothsex
  
  #random 10% dispersal philo
  dataa <- bothsex_10
  plot(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), 
       col="blue",ylim=c(0,3),xlab="Years",ylab="Mean LRS",main=paste0(c("Bothsex - Pop "),popo),pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), col="blue", lwd=2)
  
  #random 10% dispersal immi
  points(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,])[,2]~ c(12:44), col="red",pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,]), col="red", lwd=2)
  
  #females
  
  #random 10% dispersal philo
  dataa <- female_10
  plot(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), 
       col="blue",ylim=c(0,3),xlab="Years",ylab="Mean LRS",main=paste0(c("Females - Pop "),popo),pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), col="blue", lwd=2)
  
  #random 10% dispersal immi
  points(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,])[,2]~ c(12:44), col="red",pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,]), col="red", lwd=2)
  
  #males
  
  #random 10% dispersal philo
  dataa <- male_10
  plot(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), 
       col="blue",ylim=c(0,3),xlab="Years",ylab="Mean LRS",main=paste0(c("Males - Pop "),popo),pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), col="blue", lwd=2)
  
  #random 10% dispersal immi
  points(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,])[,2]~ c(12:44), col="red",pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,]), col="red", lwd=2)
  
  legend("bottomright",legend=c("10%-Philo","10%-Immi"),fill=c("blue","red"), bty='n',border=NA)
  
}
