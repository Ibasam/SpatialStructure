## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


dir_results <- "/yourfolder/results/"

#_________________ PACKAGES _________________#
library(fst)
library(imputeTS)
#devtools::install_github("ecofolio", username="seananderson")
#library(ecofolio)

#_________________ FUNCTIONS_________________#
#source("pe_mv.R")

CV <- function(x) (sd(x/mean(x)))*100

#_________________ PARAMETERS _________________#
#source("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/dataIbasam.R")
#source("/media/hdd4To/mbuoro/MetaIBASAM-Projects/FishStrat/parIbasam.R")

nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15


# #h = c(1, .95, .9, .85, .8, .75, .7) # Philopatry (homing) rates #edit al - 22/03/21 5% and 15%
# scenarioConnect=7 #scenario 1 for h=1.00, scenario 2 for h=0.95, scenario 3 pour h=0.80
# 
# # 0: control / 
# # 1: no fishing on sink populations 
# # 2: no fishing on source populations 
# scenarioFishing = 0
# 
# #frates_vec =c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
# scenarioFishingRate = 1
# 
# 
# Sc_name <- paste0("FishStrat_",scenarioConnect,scenarioFishing,scenarioFishingRate)

nParr <- list() # Nb returns
nSmolt <- list() # Nb returns
nReturns <- list() # Nb returns
Mig <- list() # Nb immigrnats
#NRet.res <- list() # Nb returns
PE <- list()
nExploitation <- list()

#LFParr <- list()
#LFSmolt <- list()
#LFReturns <- list()
Fishing_rates <- list()

# SCENARIOS
iEXPE <- args
#EXPE <- c(110, 210,220)#,120,123,124, 310, 320, 323, 324)
#EXPE <- c(120, 121, 122, 123, 124, 320, 321, 322, 323, 324)

#for (iEXPE in EXPE){ # Loop over scenario 

 
  tmp1 <- tmp2 <- tmp3 <- tmp4 <- tmp5 <- tmp6 <- tmp7 <- tmp8<-list()
  
  #_________________ DATA _________________#

  cat("Composition for ", iEXPE, "\n")
  
  pe=pe_trend=sync=sync_trend=NULL
  
  for (repmod in 1:nSIMUL){ # loop over simulations
    #cat(indice,"/",repmod,"- ")
    cat("Simulation :",repmod," / ")
    perc <- (repmod/nSIMUL)*100
    if(perc %in% (seq(0,90,10))) cat(perc,"% ~ ","\n")
    if(perc == 100) cat(perc,"%","\n")
    
    # Loading data files
    tmp0 <- NULL
    Nparr0 <- Nsmolt0 <- NRet <- N1SW <- NMSW <- Nhomers <- NIm <- Im.table <- Exploitation <- matrix(NA,nYear+nInit,npop)
    #LFparr0 <- LFsmolt0 <- LFreturns <- matrix(NA,nYear+nInit,npop)
    fishing_rates=NULL
    results=NULL
    
    
    #load(paste0(dir_results,"results_FishStrat_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".RData"))
    df=NULL
    df <- read.fst(paste0(dir_results,"results_Dispersal_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".fst"))
    #load(paste0(dir_results,"results_Dispersal_",iEXPE,"/Fishing_rates_Sc",iEXPE,"_Sim",repmod,".RData"))

        for (pop in 1:npop){ # loop over popualtions
      
      demo <- NULL
      demo <- subset(df, Pop==pop)
      #demo <- results[[pop]]$pop
      #demo <- demo[demo$year>nInit,]
      #demo$year	<- demo$year - nInit
      nyears <- max(demo$year,na.rm=TRUE)-1
      #redds <- NULL
      #redds <- results[[pop]]$redds
      
      #fishing_tmp<-NULL
      #fishing_tmp <- results[[pop]]$fishing_rate
      #fishing_rates <- rbind(fishing_rates, fishing_tmp)
      
      #-------------------------------#
      #  Composition des populations  #
      #-------------------------------#  
      
      # toto1 <- demo$year # years (from 1 to (max(demo$year)-1))
      # if(sum(is.na(toto1)==0)!=0){
      #   nReturns.tmp = NULL
      #   for (i in 1:(max(demo$year)-1)){
      #     res <- proportions.population(demo[demo$year==i,])
      #     
      #     nReturns.tmp <- c(nReturns.tmp,res$nReturns)
      # 
      #   }} else{res <- NA}
      # 
      # nReturns.table <- rbind(nReturns.table, nReturns.tmp)
      # 
      
      # Number of returns
      #toto1 <- demo$year # years (from 1 to (max(demo$year)-1))
      #toto2 <- demo$Returns # individual indicator of return ( 1, 0 otherwise)
      #toto3 <- demo$CollecID # individual indicator of return ( 1, 0 otherwise)
      #toto4 <- demo$Parr
      #toto5 <- demo$Smolt
      
      
      for (i in 1:nyears){
        
        ## PARR 0+
 #       tmp<-subset(demo,Parr==1 & date==273 & AgeRiver<1 & year==i)
 #       Nparr0[i,pop] <- nrow(tmp)
        #LFparr0[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## SMOLT 0+
 #       tmp<-subset(demo,Smolt==1 & date==90 & AgeRiver==1 & year==i)
 #       Nsmolt0[i,pop] <- nrow(tmp)
        #LFsmolt0[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        
        ## ANADROMOUS
        tmp<-subset(demo,Returns==1 & date==273 & year==i)
        NRet[i,pop] <- nrow(tmp)
        #LFreturns[i,pop] <- mean(tmp$Lf,na.rm=TRUE)
        #1SW
 #       tmp<-subset(demo,Returns==1 & AgeSea < 2 & date==273 & year==i)
 #       N1SW[i,pop] <- nrow(tmp)
        #MSW
 #       tmp<-subset(demo,Returns==1 & AgeSea >= 2 & date==273 & year==i)
 #       NMSW[i,pop] <- nrow(tmp)
        
        ## HOMERS
       tmp<-subset(demo,Returns==1 & date==273 & CollecID == pop & year==i)
       Nhomers[i,pop] <- nrow(tmp)
        
        ## IMMIGRANTS
        tmp<-subset(demo,Returns==1 & date==273 & CollecID != pop & year==i)
        NIm[i,pop] <- nrow(tmp)
        
        #Nb immigrants by population of origin
      for (pop2 in 1:npop){
         if (pop2 != pop){
           tmp<-subset(demo,Returns==1 & date==273 & CollecID == pop2 & year==i)
           Im.table[i,pop2] <- nrow(tmp)
         } else {Im.table[i,pop2] <- 0}
      }
        
        ## FISHING
 #       if (i<=nInit){
 #         Exploitation[i,pop] <- NRet[i,pop]*fishing.rates[pop,1]
 #       }else{
 #         Exploitation[i,pop] <- NRet[i,pop]*fishing.rates[pop,2]
 #       }
        
      } # end loop years
      
      tmp0[[pop]] <- Im.table
      #tmp <-subset(demo,Parr==1 & date==273 & AgeRiver<1 & year==i)
      #hist(tmp$Lf)
      
    } # end loop population
    
#    tmp1[[repmod]] <- Nparr0
#    tmp2[[repmod]] <- Nsmolt0
    
    NIm<-na_replace(NIm, 0)
    
    tmp3[[repmod]] <- list(Nhom=Nhomers, NIm=NIm, NEm=Reduce('+',tmp0), Im=tmp0)
    #tmp3[[repmod]] <- list(NIm=NIm)
    
    NRet<-na_replace(NRet, 0)
    tmp4[[repmod]] <- NRet
#    tmp5[[repmod]] <- N1SW
#    tmp6[[repmod]] <- NMSW
    
    #tmp5[[repmod]] <- LFparr0
    #tmp6[[repmod]] <- LFsmolt0
    #tmp7[[repmod]] <- LFreturns
    
#    tmp8[[repmod]] <- Exploitation
    
    
    # tmp1[[repmod]] <- NParr.table
    # tmp2[[repmod]] <- NSmolt.table
    # tmp3[[repmod]] <- list(Nhom=Nhomers.table, NIm=NIm.table, NEm=Reduce('+', tmp2), Im=tmp2)
    # tmp4[[repmod]] <- NRet.table
    #tmp6[[repmod]] <- list(Npeche=Npeche.table, Ppeche=Ppeche.table)
    
    
    ## COMPUTE PE ##
    
    #tmp <- pe_mv(NRet[11:50,], type="linear_detrended", ci = TRUE) #No initialization phase!
 #   tmp_trend <- pe_mv(NRet[11:50,], type="linear", ci = TRUE)
    
    #pe[repmod] <- tmp$pe
 #   pe_trend[repmod] <- tmp_trend$pe
    
 #   sync_trend[repmod] <- synchrony(NRet[11:50,]) #0.2635 ?a c'est bien :) car 1 synchrone, 0 asynchrone
    #sync[repmod] <- synchrony_detrended(NRet[11:50,]) #0.2635 ?a c'est bien :) car 1 synchrone, 0 asynchrone
    
    # # tmp<- pe_avg_cv(NRet.table, detrending=c("linear_detrended"))
    # tmp <- pe_mv(NRet.table, type="linear", ci = TRUE)
    # pe[repmod] <- tmp$pe
    # CV_est[repmod] <- tmp$cv_single_asset
    # CV_obs[repmod] <- tmp$cv_portfolio
    # 
    # tmp <- pe_mv(NRet.table, type = "loess_detrended", ci = TRUE)
    # petr[repmod] <- tmp$pe
    # CV_esttr[repmod] <- tmp$cv_single_asset
    # CV_obstr[repmod] <- tmp$cv_portfolio
    # # m?tapop 1.3 fois plus stable que si c'?tait une pop isol?e mais ici PAS DE CONNECTIVITE
    # # effet portfolio d? ? autre chose ici?
    #sync[repmod] <- synchrony(NRet.table) #0.2635 ?a c'est bien :) car 1 synchrone, 0 asynchrone
    
  } # end loop simul
  
  #tmp5 <- sync #cbind(pe,CV_est,CV_obs,petr,CV_esttr,CV_obstr,sync)
  #tmp5 <- cbind(pe,pe_trend,sync,sync_trend)
  tmp5 <- cbind(pe,sync)
  
#  nParr[[paste0(iEXPE)]] <- tmp1
#  nSmolt[[paste0(iEXPE)]] <- tmp2
  Mig[[paste0(iEXPE)]] <- tmp3
  nReturns[[paste0(iEXPE)]] <- tmp4
  
  #LFParr[[paste0(iEXPE)]] <- tmp5
  #LFSmolt[[paste0(iEXPE)]] <- tmp6
  #LFReturns[[paste0(iEXPE)]] <- tmp7
  
#  nExploitation[[paste0(iEXPE)]] <- tmp8
  
#  Fishing_rates[[paste0(iEXPE)]] <- fishing_rates
  
  #PE[[paste0(iEXPE)]] <- tmp5
  #Exploitation[[paste0(iEXPE)]] <- tmp6
  
  #save(nReturns=tmp4, PE=tmp5, Exploitation=tmp6, file=paste0("results/DEMOGRAPHY_",iEXPE,".RData"))
  #save(nReturns, nImm, PE, Exploitation, file="results/DEMOGRAPHY.RData")

  #} # end scenario

### Save results
#save(nReturns, Imm, PE, Exploitation, file="results/DEMOGRAPHY0.RData")
#save(nParr,nSmolt,nReturns, Mig, LFParr,LFSmolt,LFReturns,fishing.rates,nExploitation,file=paste0(dir_results,"DEMOGRAPHY",iEXPE,".RData"))
#save(nParr,nSmolt,nReturns, Mig, PE,fishing.rates,nExploitation,file=paste0(dir_results,"DEMOGRAPHY",iEXPE,".RData"))
save(nReturns,Mig,file=paste0(dir_results,"DEMOGRAPHY",iEXPE,"_50.RData"))

q('no')
