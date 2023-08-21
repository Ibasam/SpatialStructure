#########------------------------- GLOBAL DATA AND PARAMETERS -------------------------#############

nSIMUL=50 #100
npop=15
nYears=40
nInit=10
pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
dat<-dat[-1,]
load("data/random_propG_0.3.RData")

## SCENARIOS
EXPE<-c(10311111, 30311111, 50311111, 70311111, #0-30% dispersal - No diversity - simple configuration - Higher within pop SD at initialization
        10341111, 30341111, 50341111, 70341111, #0-30% dispersal - Gradual diversity - simple configuration
        10341211, 30341211, 50341211, 70341211, #0-30% dispersal - Random diversity - simple configuration
        30341221, 50341221, 70341221,           #10-30% dispersal - Random diversity - complex configuration (distance)
        30341212, 50341212, 70341212,           #10-30% dispersal - Random diversity - complex configuration (carrying capacity)
        30341222, 50341222, 70341222)           #10-30% dispersal - Random diversity - complex configuration (both)


scn_nodiv <- c(1:4)
scn_grad <- c(5:8) 
scn_rand <- c(9:12)
scn_dist <- c(13:15)
scn_size <- c(16:18)
scn_both <- c(19:21)
scn_spat_10 <- c(10,13,16,19)
scn_all_div <- c(1:12)
scn_all_spat <- c(10:21)
#optimum
optimum<-rep(0,15)

color_scn <- c("grey80","grey60","grey40","grey20",
               "darkseagreen1","mediumseagreen","forestgreen","darkolivegreen",
               "gold", "darkorange","darkorange3", "lightsalmon4",
               "red","red", "red", 
               "dodgerblue","dodgerblue","dodgerblue",
               "darkviolet", "darkviolet", "darkviolet")

color_scn_basic<-c('black',"forestgreen","orange")
color_scn_config<-c('orange',"red","dodgerblue","darkviolet")

#########--------------------------- PACKAGES AND FUNCTIONS ---------------------------#############
library(visNetwork)
library(reshape2)
library(imputeTS)
library(ggplot2)
cv <- function(x){
  sd(x)/abs(mean(x))
}
se <- function(x) sd(x,na.rm=TRUE) / sqrt(length(x)) # stabdard error
#The standard error is the standard deviation of the mean in repeated samples from a population.

#########--------------------------- LOAD RESULTS FROM SIMULATIONS ---------------------------#############

res_1SW_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/PHENOGENOTYPE",EXPE[iEXPE],"_50.RData"))
  res_1SW_scn[[iEXPE]]<-res_1SW
}

nReturns_scn <- list()
Mig_scn <- list()
for (iEXPE in 1:length(EXPE)) {
  load(paste0("results/DEMOGRAPHY",EXPE[iEXPE],"_50.RData"))
  nReturns_scn[[iEXPE]]<-nReturns
  Mig_scn[[iEXPE]]<-Mig
}


#########---------------------------- CODE FOR FIGURES ----------------------------#############

##########################################################################################
########## COMPUTE METRICS OF LOCAL TRAIT MISMATCH AND EVOLUTIONARY RATES ################
##########################################################################################

ltm_vit<-array(,dim=c(npop,length(EXPE)))
ltm_ratio<-array(,dim=c(npop,length(EXPE)))
ltm_50<-array(,dim=c(npop,length(EXPE)))
ltm_40<-array(,dim=c(npop,length(EXPE)))
ltm_30<-array(,dim=c(npop,length(EXPE)))
ltm_init<-array(,dim=c(npop,length(EXPE)))

cv_years<-array(,dim=c(nInit+nYears,npop,length(EXPE)))
cv_years_immi<-array(,dim=c(nInit+nYears,npop,length(EXPE)))

for (pop in 1:npop) { 
  geno_1SW_tot<-list()
  for (scn in c(1:length(EXPE))){
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  for (scn in 1:length(EXPE)) {
    mean_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), mean)[1:(nInit+nYears),2]
    cv_gene<-aggregate(exp(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop]), by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), cv)[1:(nInit+nYears),2]

    ltm_init[pop,scn]<-mean(mean_gene[1:5], na.rm=T)
    ltm_50[pop,scn]<-mean(mean_gene[(nInit+nYears-5):(nInit+nYears)], na.rm=T)
    ltm_40[pop,scn]<-mean(mean_gene[(nInit+nYears-15):(nInit+nYears-10)], na.rm=T)
    ltm_30[pop,scn]<-mean(mean_gene[(nInit+nYears-25):(nInit+nYears-20)], na.rm=T)
    
    if(ltm_50[pop,scn] < ltm_init[pop,scn]){
      ltm_vit[pop,scn]<- min(which(abs(mean_gene)[1:(nInit+nYears)] <= (mean(abs(mean_gene)[(nInit+nYears-5):(nInit+nYears)], na.rm=T))))#+mean(abs(mean_gene)[(nInit+nYears-5):(nInit+nYears)], na.rm=T)*0.01)))
    } else {
      ltm_vit[pop,scn]<- min(which(abs(mean_gene)[1:(nInit+nYears)] >= (mean(abs(mean_gene)[(nInit+nYears-5):(nInit+nYears)], na.rm=T))))#+mean(abs(mean_gene)[(nInit+nYears-5):(nInit+nYears)], na.rm=T)*0.01)))
    }
    ltm_ratio[pop,scn]<-abs(mean(mean_gene[1:5], na.rm=T)-mean(mean_gene[(nInit+nYears-5):(nInit+nYears)], na.rm=T))
    ltm_ratio[pop,scn]<-ltm_ratio[pop,scn]/ltm_vit[pop,scn]

    cv_years[,pop,scn]<-cv_gene

    if (scn %in% c(2:4,6:8,10:21)) {
      cv_gene_immi<-aggregate(exp(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID!=pop]), by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), cv)[1:(nInit+nYears),2]
      cv_years_immi[,pop,scn]<-cv_gene_immi
    }
  }
}


##########################################################################################
######################## COMPUTE METRICS OF POPULATIONS DEMOGRAPHY #######################
##########################################################################################

#Population size
Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(EXPE),npop))
for (pop in 1:npop){
  for (scn in c(1:length(EXPE))){
    for (simul in 1:nSIMUL){
      Npops[,simul,scn,pop] <- nReturns_scn[[scn]][[1]][[simul]][1:(nYears+nInit),pop] #nb returns metapop per year, simu and scenario
    }
  }
}

#Proportion of Immigrants
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

##########################################################################################
######################## COMPUTE METRICS OF METAPOPULATION DEMOGRAPHY ####################
##########################################################################################

nMetapop_size_init <- nMetapop_size_5 <- nMetapop_size_relative <-nMetapop_cv<-nMetapop_var <-nMetapop_mean<- array(,dim=c(nSIMUL,length(EXPE)))
Npops_mean<-Npops_var<-array(,dim=c(npop,nSIMUL,length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    nMetapop_size_init[simul,scn]<-median(rowSums(nReturns_scn[[scn]][[1]][[simul]])[1:2])
    nMetapop_size_5[simul,scn]<-median(rowSums(nReturns_scn[[scn]][[1]][[simul]])[45:50])
    nMetapop_size_relative[simul,scn]<-median(rowSums(nReturns_scn[[scn]][[1]][[simul]])[45:50])/median(rowSums(nReturns_scn[[scn]][[1]][[simul]])[1:2])
    nMetapop_cv[simul,scn]<-cv(rowSums(nReturns_scn[[scn]][[1]][[simul]])[1:50])
    nMetapop_var[simul,scn]<-var(rowSums(nReturns_scn[[scn]][[1]][[simul]])[1:50])
    nMetapop_mean[simul,scn]<-mean(rowSums(nReturns_scn[[scn]][[1]][[simul]])[1:50])
  }
}


################################################################
### BOXPLOT METRICS METAPOPULATION SCALE (Fig. 3, Fig. S12)  ###
################################################################

######## EVOLUTION #######

metapop_evoPlot <- function(metric, ylim, ylab, scenarios){
  
  if(scenarios == "diversity") {
    metapop_scenarios <- list(metric[,scn_nodiv], metric[,scn_grad], metric[,scn_rand])
    gap <- c(-.15,0,.15)
    plot(NULL, xlim=c(0.5,4.5), ylim=ylim, xlab="Dispersal rate",ylab=ylab, xaxt='n', cex.axis=1.2, cex.lab=1.2)
    mtext(c("0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4, cex=1.2)
    
    for (i in 1:length(metapop_scenarios)){
      points(rep(c(1:4+gap[i]),each=15),metapop_scenarios[[i]],pch=4, cex=0.5,col=color_scn_basic[i])
      points(c(1:4+gap[i]),apply(metapop_scenarios[[i]],2,median),pch=19, cex=1.5,col=color_scn_basic[i])
      segments(c(1:4+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.25),
               c(1:4+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.75),lwd=3,col=color_scn_basic[i])
      
    }
    legend("topright", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=color_scn_basic, bty="n", ncol=1)
  }else{
    gap<-c(-.3,-.1,.1,.3)
    metapop_scenarios <- list(metric[,scn_rand[-1]], metric[,scn_dist], metric[,scn_size], metric[,scn_both])
    
    plot(NULL, xlim=c(0.5,3.5), ylim=ylim, xlab="Dispersal rate",ylab=ylab, xaxt='n', cex.axis=1.2, cex.lab=1.2)
    mtext(c("10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:3, cex=1.2)
    
    for (i in 1:length(metapop_scenarios)){
      points(rep(c(1:3+gap[i]),each=15),metapop_scenarios[[i]],pch=4, cex=0.5,col=color_scn_config[i])
      points(c(1:3+gap[i]),apply(metapop_scenarios[[i]],2,median),pch=19, cex=1.5,col=color_scn_config[i])
      segments(c(1:3+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.25),
               c(1:3+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.75),lwd=3,col=color_scn_config[i])
    }
    legend("topright", cex=.8, legend=c("Simple config.","Distances","Carr. capacity", "Dist. + Carr. cap."), col=2:3, border=NA,fill=color_scn_config, bty="n", ncol=1)
  }
}

##Local trait mismatch
y=50
metapop_evoPlot(eval(parse(text=paste0("ltm_","",y))), c(0,0.15),"Final local trait mismatch", "diversity")
#mtext(paste0(y,"yrs"), line=2)
abline(h=0.15,lty=2,col="red")

y=50
metapop_evoPlot(eval(parse(text=paste0("ltm_","",y))), c(0,0.15),"Final local trait mismatch", "config")
#mtext(paste0(y,"yrs"), line=2)
abline(h=0.15,lty=2,col="red")

##Evolutionary rate
metapop_evoPlot(ltm_ratio*1000, c(0,8),"Evolutionary rate", "diversity")

metapop_evoPlot(ltm_ratio*1000, c(0,8),"Evolutionary rate", "config")


###### DEMOGRAPHY ######

metapop_demoPlot <- function(metric, ylim, ylab, scenarios){
  
  if(scenarios == "diversity") {
    metapop_scenarios <- list(metric[,scn_nodiv], metric[,scn_grad], metric[,scn_rand])
    gap <- c(-.15,0,.15)
    plot(NULL, xlim=c(0.5,4.5), ylim=ylim, xlab="Dispersal rate",ylab=ylab, xaxt='n', cex.axis=1.2, cex.lab=1.2)
    mtext(c("0%","10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:4, cex=1.2)
    
    for (i in 1:length(metapop_scenarios)){
      points(c(1:4+gap[i]),apply(metapop_scenarios[[i]],2,median),pch=19, cex=1.5,col=color_scn_basic[i])
      segments(c(1:4+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.25),
               c(1:4+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.75),lwd=3,col=color_scn_basic[i])
    }
    legend("topright", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=color_scn_basic, bty="n", ncol=1)

  }else{
    gap<-c(-.3,-.1,.1,.3)
    metapop_scenarios <- list(metric[,scn_rand[-1]], metric[,scn_dist], metric[,scn_size], metric[,scn_both])
    
    plot(NULL, xlim=c(0.5,3.5), ylim=ylim, xlab="Dispersal rate",ylab=ylab, xaxt='n', cex.axis=1.2, cex.lab=1.2)
    mtext(c("10%","20%","30%"), side = 1, line = 1, outer = FALSE, at = 1:3, cex=1.2)
    
    for (i in 1:length(metapop_scenarios)){
      points(c(1:3+gap[i]),apply(metapop_scenarios[[i]],2,median),pch=19, cex=1.5,col=color_scn_config[i])
      segments(c(1:3+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.25),
               c(1:3+gap[i]),apply(metapop_scenarios[[i]],2,quantile,probs=.75),lwd=3,col=color_scn_config[i])
    }
    legend("topright", cex=.8, legend=c("Simple config.","Distances","Carr. capacity", "Dist. + Carr. cap."), col=2:3, border=NA,fill=color_scn_config, bty="n", ncol=1)
  }
}

#Metapopulation size
metapop_demoPlot(nMetapop_size_5, c(500,2000),"Metapopulation size", "diversity")
#metapop_demoPlot(nMetapop_size_relative, c(0.2,1),"Relative metapopulation size", "diversity")

metapop_demoPlot(nMetapop_size_5, c(500,2000),"Metapopulation size", "config")

#Metapopulation CV
metapop_demoPlot(nMetapop_cv, c(0,0.3),"Metapopulation CV", "diversity")

metapop_demoPlot(nMetapop_cv, c(0,0.3),"Metapopulation CV", "config")




################################################################################
###### TEMPORAL EVOLUTION OF METRICS BY POPULATION SCN DIVERSITY (Fig. 4) ######
################################################################################

######## GROWTH POTENTIAL #######

#separated, by gG init rather than pop number
scenarios <- list(scn_nodiv, scn_grad, scn_rand)
gG_init<-round(seq(from=0.3, to=0, length.out=npop),2)
for(j in 1:length(scenarios)) {
  par(mfrow=c(3,3))
  for (i in 1:npop) {
    plot(NULL, xlim=c(0,50), ylim=c(-.01,.27), xlab="Years",ylab="Growth potential", main=gG_init[i], cex.lab=1.2,cex.axis=1.2)
    abline(h=0, lty=2)
    # if(j=1 && i!=8){
    #   
    # } else {
      for (scn in scenarios[[j]]) { 
        if(scn>8){ #scn random
          pop<-rev(order(random_propG))[i]
        }else{
          pop<-i
        }
        geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
        for (simul in 2:nSIMUL) {
          geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
        }
        
        median_gene<-aggregate(geno_1SW$gG[geno_1SW$CollecID==pop], by=list(geno_1SW$years[geno_1SW$CollecID==pop]), median)[1:50,2]
        lw<-loess(median_gene~c(1:50),span=0.2)
        lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
        if(scn %in% c(2:4,6:8,10:12) ) {
          #if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
          median_gene_immi<-aggregate(geno_1SW$gG[geno_1SW$CollecID!=pop], by=list(geno_1SW$years[geno_1SW$CollecID!=pop]), median)[1:50,2]
          lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
          if(scn %in% c(2:4)) {
            lines((1:49), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[1:49,pop,scn]),lty=3)
          }else{
            lines((1:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[1:50,pop,scn]),lty=3)
          }
        }
      }
    #}
  }
}

  
######## POPULATION DEMOGRAPHY #######
scenarios <- list(scn_nodiv, scn_grad, scn_rand)

par(mfrow=c(3,3))

for (i in 1:npop) {
  plot(NULL, xlim=c(0,(nYears+nInit)), ylim=c(0,1), xlab="Years",ylab="Relative abundance",main=gG_init[i], cex.lab=1.2,cex.axis=1.2) #500
  
  for(j in 1:length(scenarios)) {
    if(j==1 && i!=8) {

    }else{
      for (scn in scenarios[[j]]) { 
        if(scn>8){
          pop<-rev(order(random_propG))[i]
        }else{
          pop<-i
        }
        median_npops<- apply(Npops[,,scn,pop],1,median)
        lw<-loess(median_npops~c(1:(nYears+nInit)),span=0.2)
        lines((1:(nYears+nInit)), lw$fitted/lw$fitted[1], pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
      }
    }
  }
}  


######################################################################################
### TEMPORAL EVOLUTION OF METRICS BY POPULATION SCN SPATIAL CONFIGURATION (Fig. 5) ###
######################################################################################
gG_random_init<-round(random_propG,2)

######## GROWTH POTENTIAL #######

par(mfrow=c(3,3))
for (pop in 1:npop) {
  geno_1SW_tot<-list()
  for (scn in scn_spat_10){ #define scenarios
    geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
    for (simul in 2:nSIMUL) {
      geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
    }
    geno_1SW_tot[[scn]] <- geno_1SW
  }
  
  plot(NULL, xlim=c(0,(nInit+nYears)), ylim=c(-.01,.27), xlab="Years",ylab="Genotypic Growth potential", main=paste0("Pop ",pop," - LTM ",gG_random_init[pop]), cex.lab=1.2,cex.axis=1.2)
  abline(h=0, lty=2)
  
  for (scn in scn_spat_10) { #define scenarios here too
    median_gene<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID==pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID==pop]), median)[1:(nInit+nYears),2]
    lw<-loess(median_gene~c(1:(nInit+nYears)),span=0.2)
    lines((1:(nInit+nYears)), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
    
    #if(scn %in% c(2:4,6:8,10:12) ) {
    if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
      median_gene_immi<-aggregate(geno_1SW_tot[[scn]]$gG[geno_1SW_tot[[scn]]$CollecID!=pop], by=list(geno_1SW_tot[[scn]]$years[geno_1SW_tot[[scn]]$CollecID!=pop]), median)[1:(nInit+nYears),2]
      lw_immi<-loess(median_gene_immi~c(1:(nInit+nYears)),span=0.2)
      lines((1:(nInit+nYears)), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[1:(nInit+nYears),pop,scn]),lty=3)
    }
  }
  #legend("topright", cex=.8, legend=c("Simple config.","Distances","Carr. capacity", "Dist. + Carr. cap."), col=2:3, border=NA,fill=color_scn_config, bty="n", ncol=1)
}


######## POPULATION DEMOGRAPHY #######

par(mfrow=c(3,3))
for (pop in 1:npop){
  Npops <- array(,dim=c((nYears+nInit),nSIMUL,length(EXPE)))
  for (scn in 1:length(EXPE)){
    for (simul in 1:nSIMUL){
      Npops[,simul,scn] <- nReturns_scn[[scn]][[1]][[simul]][1:(nYears+nInit),pop] #nb returns metapop per year, simu and scenario
    }
  }
  
  plot(NULL, xlim=c(0,(nYears+nInit)), ylim=c(0,1.5), xlab="Years",ylab="Relative abundance",main=paste0("Pop ",pop," - LTM ",gG_random_init[pop]),, cex.lab=1.2,cex.axis=1.2)

  #lw_init<-NULL
  for (scn in scn_spat_10) { #choice of scenarios to plot
    median_npops<- apply(Npops[,,scn],1,median)
    lw<-loess(median_npops~c(1:(nYears+nInit)),span=0.2)
    #lw_init<-c(lw_init,lw$fitted[1])
    lines((1:(nYears+nInit)), lw$fitted/lw$fitted[1], pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
  }
  #legend("topright", cex=.8, legend=c("Simple config.","Distances","Carr. capacity", "Dist. + Carr. cap."), col=2:3, border=NA,fill=color_scn_config, bty="n", ncol=1)
}






################################################# SUPPLEMENTARY MATERIALS ##################################################



###############################################################
####### Fig.S2 RATIO IMMIGRANTS/EMIGRANTS BY POPULATION  ######
###############################################################

I <- P.im <- array(,dim=c(nSIMUL,npop, length(EXPE)))
ratio2 <- array(,dim=c((nYears+nInit),npop,nSIMUL, length(EXPE)))
for (scn in c(scn_nodiv[-1],scn_grad[-1],scn_all_spat)){
  for (simul in 1:nSIMUL){
    ratio<-Mig_scn[[scn]][[1]][[simul]]$NIm/Mig_scn[[scn]][[1]][[simul]]$NEm
    ratio[which(ratio=="Inf")]<-NA
    I[simul, ,scn] <- log(apply(ratio[45:50,],2,median, na.rm=T)) #5last years #mean , na.rm=T
  }
}

gap<-c(-0.3,-0.15,  0.0,  0.15, -0.3, 0,  0.0,  0.15, -0.3, 0.15,  0.0,  0.15, 0,  0.0,  0.15, 0.15,  0.0,  0.15, 0.3,  0.0,  0.15)
plot(NULL, xlim=c(0,(npop+1)), ylim=c(-1,2), xlab="Population",ylab="Ratio Immigrants/Emigrants (log)", xaxt='n')
mtext(1:npop, side = 1, line = 1, outer = FALSE, at = 1:15)
abline(h=0,col="black",lty=2)
color_scn_here=c("","black","","","","forestgreen","","","","orange")
for (scn in c(scn_nodiv[2],scn_grad[2],scn_rand[2])){ ##10% dispersal
  for (pop in 1:npop) {
    points(pop+gap[scn], median(I[,pop,scn],na.rm=TRUE), pch=20, cex=1.8,col=color_scn_here[scn]);
    segments(pop+gap[scn], quantile(I[,pop,scn],probs=.25,na.rm=TRUE),pop+gap[scn],quantile(I[,pop,scn],probs=.75,na.rm=TRUE),col=color_scn_here[scn], lwd=2)
  }
}

gap<-c(-0.3,-0.15,  0.0,  0.15, -0.3, 0,  0.0,  0.15, -0.3, -0.15,  0.0,  0.15, 0,  0.0,  0.15, 0.15,  0.15,  0.15, 0.3,  0.0,  0.15)
plot(NULL, xlim=c(0,(npop+1)), ylim=c(-1,2), xlab="Population",ylab="Ratio Immigrants/Emigrants (log)", xaxt='n')
mtext(1:npop, side = 1, line = 1, outer = FALSE, at = 1:15)
abline(h=0,col="black",lty=2)
for (scn in scn_spat_10){ ##10% dispersal
  for (pop in 1:npop) {
    points(pop+gap[scn], median(I[,pop,scn],na.rm=TRUE), pch=20, cex=1.8,col=color_scn[scn]);
    segments(pop+gap[scn], quantile(I[,pop,scn],probs=.25,na.rm=TRUE),pop+gap[scn],quantile(I[,pop,scn],probs=.75,na.rm=TRUE),col=color_scn[scn], lwd=2)
  }
}


################################
### Fig.S3 NETWORKS ALL SCN ####
################################

###########  Migrants flows

a <- array(,dim=c(npop,npop,nSIMUL, length(EXPE)))
for (scn in c(scn_nodiv[-1],scn_grad[-1],scn_all_spat)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      a[pop,,simul,scn]<-apply(Mig_scn[[scn]][[1]][[simul]]$Im[[pop]][1:5,],2,mean)
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
  arr[[scn]] <- array( unlist(a[,,,scn]) , c(npop,npop,100) )
  mflow[[scn]]<-apply( arr[[scn]] , 1:2 , mean )#mean simulations
  colnames(mflow[[scn]])<-pops
  rownames(mflow[[scn]])<-pops
}

########### Populations size
Parpop<-array(dim=c(nSIMUL, npop, length(EXPE)))
Parpop2<-array(dim=c(npop, length(EXPE)))
for (scn in 1:length(EXPE)) {
  for (simul in 1:nSIMUL) {
    for (pop in 1:npop) {
      Parpop[simul, ,scn] <- colMeans(nReturns_scn[[scn]][[1]][[simul]][1:5,], na.rm=T) #5last years
    }
  }
  Parpop2[,scn]<-colMeans(Parpop[,,scn])
}

visNetworkPlot <- function(scenario){
  
  mflow2<-melt(mflow[[scenario]]) #for absolute values 
  mflow3<-mflow2[c(2,1,3)] #for absolute values
  
  mflow3[,1]<-as.factor(mflow3[,1])
  mflow3[,2]<-as.factor(mflow3[,2])
  
  nodes<-data.frame(id=c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15"),
                    pop=pops,
                    type=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"), #control
                    size=Parpop2[,scenario]
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
  vis.nodes$size   <- nodes$size/5#5 # Node size
  vis.nodes$borderWidth <- 2 # Node border width
  
  colors<-c('#94003a', '#9e2045', '#a83351', '#b2445e', '#bc546b', '#c56478', '#ce7385', '#d78293', '#df91a0', '#e7a1af', '#eeb1bd', '#f4c0cc', '#f9d1db', '#fce2ea', '#fcf3fa')
  
  if(scenario %in% scn_nodiv) {
    #for scn no diversity
    gradient_color <- rep(colors[8],npop) 
  }
  if(scenario %in% scn_grad) {
    #for scn gradual diversity
    gradient_color <- colors
  }
  if(scenario %in% scn_spat_10) {
    #for scn random diversity
    growth<-c(0.23571429, 0.21428571, 0.00000000, 0.08571429, 0.10714286, 0.25714286, 
              0.04285714, 0.17142857, 0.12857143, 0.15000000, 0.06428571, 0.19285714, 
              0.02142857, 0.30000000, 0.27857143)
    gradient_color<-rev(colors)[rank(growth)]
  }
  
  vis.nodes$color.background <- gradient_color#colors #rep(colors[8],15) #gradient_color_random #rev(gradient_color_func(15))[as.numeric(as.character(nodes$type))]
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
  if(scenario %in% c(1:12,16:18)) {
    lat2<-seq(max(lat),min(lat),by=-(max(lat)-min(lat))/14)
    lon2<-c(seq(max(lon),min(lon),by=-(max(lon)-min(lon))/7),seq(min(lon),max(lon),by=(max(lon)-min(lon))/7))
    lon2<-lon2[-9]
    lat<-lat2
    lon<-lon2
    #plot(lon2, lat2)
  }
  
  #here put coordinates wanted
  vis.nodes$x<- lon*1000
  vis.nodes$y <- -lat*1000
  
  a<-visNetwork(vis.nodes, vis.links)
  a<-visEdges(a, arrows=list(to=list(enable=T, scaleFactor=1.5)),color = list(color = "black", highlight = "red"), smooth = list(enabled = TRUE, type = "diagonalCross"))
  a<-visNodes(a,fixed = TRUE,physics=T, font=list(color="black", size=0))
  a
  
}

visNetworkPlot(2)
visNetworkPlot(6)
visNetworkPlot(10)
visNetworkPlot(13)
visNetworkPlot(16)
visNetworkPlot(19)

##################################################################
########## TEMPORAL EVOLUTION OF TRAITS (Fig. S4-S10) ############
##################################################################


gG_init<-round(seq(from=0.3, to=0, length.out=npop),2)

traitPlot <- function(trait,ylim,scenarios){
  par(mfrow=c(3,3))
  if(scenarios == "diversity") {
    scenarios <- list(scn_nodiv, scn_grad, scn_rand)
    
    for (i in 1:npop) {
      plot(NULL, xlim=c(0,50), ylim=ylim, xlab="Years",ylab=trait, main=gG_init[i], cex.lab=1.2,cex.axis=1.2)
      abline(h=0, lty=2)
      
      for(j in 1:length(scenarios)) {
        if(j==1 && i!=8) {
        }else{
          for (scn in scenarios[[j]]) { #choice of scenarios to plot
            if(scn>8){
              pop<-rev(order(random_propG))[i]
            }else{
              pop<-i
            }
            geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
            for (simul in 2:nSIMUL) {
              geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
            }
            
            median_gene<-aggregate(geno_1SW[geno_1SW$CollecID==pop,trait], by=list(geno_1SW$years[geno_1SW$CollecID==pop]), median)[1:50,2]
            lw<-loess(median_gene~c(1:50),span=0.2)
            lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
            if(scn %in% c(2:4,6:8,10:12) ) {
              #if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
              median_gene_immi<-aggregate(geno_1SW[geno_1SW$CollecID!=pop,trait], by=list(geno_1SW$years[geno_1SW$CollecID!=pop]), median)[1:50,2]
              lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
              if(scn %in% c(2:4)) {
                lines((1:49), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[1:49,pop,scn]),lty=3)
              }else{
                lines((1:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[1:50,pop,scn]),lty=3)
              }
            }
          }
        }
      } 
    }
  }else{
    for (pop in 1:npop) {
      plot(NULL, xlim=c(0,50), ylim=ylim, xlab="Years",ylab=trait, main=pop, cex.lab=1.2,cex.axis=1.2)
      abline(h=0, lty=2)
      
      for (scn in scn_spat_10){ #define scenarios
        geno_1SW <- res_1SW_scn[[scn]][[1]][[pop]]
        for (simul in 2:nSIMUL) {
          geno_1SW <- rbind(geno_1SW, res_1SW_scn[[scn]][[simul]][[pop]])
        }
        
        median_gene<-aggregate(geno_1SW[geno_1SW$CollecID==pop,trait], by=list(geno_1SW$years[geno_1SW$CollecID==pop]), median)[1:50,2]
        lw<-loess(median_gene~c(1:50),span=0.2)
        lines((1:50), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=1.5)
        
        if(scn >9 ) { #plot immigrants (for scn with dispersal rates > 0)
          median_gene_immi<-aggregate(geno_1SW[geno_1SW$CollecID!=pop,trait], by=list(geno_1SW$years[geno_1SW$CollecID!=pop]), median)[1:50,2]
          lw_immi<-loess(median_gene_immi~c(1:50),span=0.2)
          lines((1:50), lw_immi$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=10*median(P.im.median[1:50,pop,scn]),lty=3)
        }
      }
    }
  }
}
traitPlot("gG",c(-.01,.27),"diversity")
traitPlot("pG",c(-.01,.27),"diversity")
traitPlot("gFmid1",c(1,2),"diversity")
traitPlot("gFmid3",c(35,55),"diversity")
traitPlot("gFmid4",c(70,130),"diversity")

traitPlot("gG",c(-.01,.27),"config")
traitPlot("pG",c(-.01,.27),"config")
traitPlot("gFmid1",c(1,2),"config")
traitPlot("gFmid3",c(35,55),"config")
traitPlot("gFmid4",c(70,130),"config")




###############################################################################
####### RELATIONSHIPS BETWEEN METRICS AND IMMIGRANTS FEATURES (Fig. S11) ######
###############################################################################

## COMPUTE METAPOPULATION TRAIT MISMATCH ##

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
dev.off()
ltm_immi <- array(,dim=c(npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  ltm_immi[,scn] <- apply(diff_trait2[45:50,,scn],2,mean)
}
data_test <- cbind(melt(ltm_50),melt(ltm_immi))
colnames(data_test) <- c("pop","scn","ltm","pop2","scn2","ltm_immi")
data_test$scn_type <- NA
for (i in 1:nrow(data_test)) {
  if (data_test$scn[i] %in% c(2:4)) {
    data_test$scn_type[i] <- "1.No_diversity"
  }
  if (data_test$scn[i] %in% c(6:8)) {
    data_test$scn_type[i] <- "2.Gradual"
  }
  if (data_test$scn[i] %in% c(10:12)) {
    data_test$scn_type[i] <- "3.Random"
  }
  if (data_test$scn[i] %in% c(13:15)) {
    data_test$scn_type[i] <- "4.Distances"
  }
  if (data_test$scn[i] %in% c(16:18)) {
    data_test$scn_type[i] <- "5.Carr_capacity"
  }
  if (data_test$scn[i] %in% c(19:21)) {
    data_test$scn_type[i] <- "6.Dist+Carr_Cap."
  }
}

ggplot(data_test, aes(x = ltm_immi, y = ltm, color = scn_type)) +
  geom_point() +
  stat_ellipse(geom = "polygon",
               aes(fill = scn_type), 
               alpha = 0.1) +
  theme_bw() +
  scale_fill_manual(values=c("black","forestgreen","orange","red","dodgerblue","darkviolet"))+
  scale_color_manual(values=c("black","forestgreen","orange","red","dodgerblue","darkviolet"))+
  xlab("Metapopulation trait mismatch") +
  ylab("Local trait mismatch") +
  theme(axis.text = element_text(size=11), axis.title= element_text(size=13))


#Evolutionary rate function of proportion of immigrants
##with ellipses
prop_immi <- array(,dim=c(npop,length(EXPE)))
for (scn in 1:length(EXPE)) {
  prop_immi[,scn] <- apply(P.im.median[1:50,,scn],2,mean)
}
data_test <- cbind(melt(ltm_ratio),melt(prop_immi))
colnames(data_test) <- c("pop","scn","ltm_ratio","pop2","scn2","prop_immi")
data_test$ltm_ratio<-data_test$ltm_ratio*1000
#data_test<-subset(data_test, data_test$pop !=6)
#data_test<-subset(data_test, data_test$pop !=7)
#data_test<-subset(data_test, data_test$pop !=8)
data_test$scn_type <- NA
for (i in 1:nrow(data_test)) {
  if (data_test$scn[i] %in% c(2:4)) {
    data_test$scn_type[i] <- "1.No_diversity"
  }
  if (data_test$scn[i] %in% c(6:8)) {
    data_test$scn_type[i] <- "2.Gradual"
  }
  if (data_test$scn[i] %in% c(10:12)) {
    data_test$scn_type[i] <- "3.Random"
  }
  if (data_test$scn[i] %in% c(13:15)) {
    data_test$scn_type[i] <- "4.Distances"
  }
  if (data_test$scn[i] %in% c(16:18)) {
    data_test$scn_type[i] <- "5.Carr_capacity"
  }
  if (data_test$scn[i] %in% c(19:21)) {
    data_test$scn_type[i] <- "6.Dist+Carr_Cap."
  }
}

try<- subset(data_test, !is.na(data_test$scn_type))
ggplot(try, aes(x = prop_immi, y = ltm_ratio, color = scn_type)) +
  geom_point() +
  stat_ellipse(geom = "polygon",
               aes(fill = scn_type), 
               alpha = 0.1) +
  theme_bw() +
  scale_fill_manual(values=c("black","forestgreen","orange","red","dodgerblue","darkviolet"))+
  scale_color_manual(values=c("black","forestgreen","orange","red","dodgerblue","darkviolet"))+
  xlab("Proportion of immigrants") +
  ylab("Evolutionary rate") +
  theme(axis.text = element_text(size=11), axis.title= element_text(size=13))





#################################################################
######### TEMPORAL AVERAGE CV GROWTH POTENTIAL (Fig. S13) #########
#################################################################

dev.off()
metapop_CVPlot <- function(metric, ylim, ylab, scenarios, immi) {
  if(scenarios == "diversity") {
    scenarios <- list(scn_nodiv, scn_grad, scn_rand)
    if(immi == TRUE) {
      plot(NULL, xlim=c(0,(nInit+nYears)), ylim=ylim, xlab="Years",ylab=ylab, cex.lab=1.2,cex.axis=1.2)
      scenarios <- list(scn_nodiv[-1], scn_grad[-1], scn_rand[-1])
    }else{
      plot(NULL, xlim=c(0,(nInit+nYears)), ylim=ylim, xlab="Years",ylab=ylab, main="", cex.lab=1.2,cex.axis=1.2)
    }
    legend("topright", cex=.8, legend=c("No diversity","Gradual diversity","Random diversity"), col=2:3, border=NA,fill=color_scn_basic, bty="n", ncol=1)
    
  }else{
    scenarios <- list(scn_rand[-1][2], scn_dist[2], scn_size[2], scn_both[2])
    if(immi == TRUE) {
      plot(NULL, xlim=c(0,(nInit+nYears)), ylim=ylim, xlab="Years",ylab=ylab, main="20%disp", cex.lab=1.2,cex.axis=1.2)
    }else{
      plot(NULL, xlim=c(0,(nInit+nYears)), ylim=ylim, xlab="Years",ylab=ylab, main="20%disp", cex.lab=1.2,cex.axis=1.2)
    }
    legend("topright", cex=.8, legend=c("Simple config.","Distances","Carr. capacity", "Dist. + Carr. cap."), col=2:3, border=NA,fill=color_scn_config, bty="n", ncol=1)
    
  }
  for (j in 1:length(scenarios)) {
    for (scn in scenarios[[j]]) { #define scenarios here too
      lw<-loess(apply(metric[,,scn],1,median)~c(1:(nInit+nYears)),span=0.3)
      if(scn %in% c(2:4) & immi==TRUE) {
        lines((1:(nInit+nYears-1)), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=3)
      }else {
        lines((1:(nInit+nYears)), lw$fitted, pch=19,cex=1.5, col=color_scn[scn], lwd=3)
      }
    }
  }
}

metapop_CVPlot(cv_years, c(0.06,.1),"Averaged local CV growth potential", "diversity",FALSE)
metapop_CVPlot(cv_years, c(0.06,.1),"Averaged local CV growth potential", "config",FALSE)
metapop_CVPlot(cv_years_immi, c(0.065,.11),"Averaged immigrants CV growth potential", "diversity",TRUE)
metapop_CVPlot(cv_years_immi, c(0.065,.11),"Averaged immigrants CV growth potential", "config",TRUE)




####################################
### Fig.S14 REPRODUCTIVE SUCCESS ###
####################################

##focus on scn dispersal 10% random diversity simple spatial configuration

load(paste0("results/FITNESS",30341211,"_50.RData"))

female_simul<-male_simul<-bothsex<-NULL
selection_intensity_fem<-array(,dim=c(45,npop,nSIMUL))
selection_intensity_mal<-array(,dim=c(45,npop,nSIMUL))
selection_intensity_bothsex<-array(,dim=c(45,npop,nSIMUL))

for (simul in 1:nSIMUL) {
  cat("Simulation :",simul," / ")
  test<-NULL
  for (pop in 1:npop) {
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
for (popo in 1:npop) {
  
  #bothsex
  
  #random 10% dispersal philo
  dataa <- bothsex_10
  plot(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), 
       col="blue",ylim=c(0,3),xlab="Years",ylab="Mean LRS",main=paste0(c("Bothsex - Pop "),popo),pch=20, cex.axis=1.2, cex.lab=1.2)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), col="blue", lwd=2)
  
  #random 10% dispersal immi
  points(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,])[,2]~ c(2:44), col="red",pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,]), col="red", lwd=2)
  
  #females
  
  #random 10% dispersal philo
  dataa <- female_10
  plot(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), 
       col="blue",ylim=c(0,3),xlab="Years",ylab="Mean LRS",main=paste0(c("Females - Pop "),popo),pch=20, cex.axis=1.2, cex.lab=1.2)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), col="blue", lwd=2)
  
  #random 10% dispersal immi
  points(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,])[,2]~ c(2:44), col="red",pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,]), col="red", lwd=2)
  
  #males
  
  #random 10% dispersal philo
  dataa <- male_10
  plot(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), 
       col="blue",ylim=c(0,3),xlab="Years",ylab="Mean LRS",main=paste0(c("Males - Pop "),popo),pch=20, cex.axis=1.2, cex.lab=1.2)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID==popo & dataa$years<45,]), col="blue", lwd=2)
  
  #random 10% dispersal immi
  points(aggregate(fitness ~years, mean, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,])[,2]~ c(2:44), col="red",pch=20)
  abline(lm(fitness ~years, data=dataa[dataa$pop==popo & dataa$CollecID!=popo & dataa$years<45,]), col="red", lwd=2)
  
  legend("bottomright",legend=c("Philopatric","Immigrants"),fill=c("blue","red"), bty='n',border=NA)
  
}
