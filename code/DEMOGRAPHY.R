## remove (almost) everything in the working environment.
rm(list = ls())
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  args <- 1
}


dir_results <- "results/"

#_________________ PACKAGES _________________#
library(fst)
library(imputeTS)

#_________________ FUNCTIONS_________________#
CV <- function(x) (sd(x/mean(x)))*100

#_________________ PARAMETERS _________________#

nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15


nParr <- list() # Nb returns
nSmolt <- list() # Nb returns
nReturns <- list() # Nb returns
Mig <- list() # Nb immigrnats
PE <- list()
nExploitation <- list()


# SCENARIOS
iEXPE <- args

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
    results=NULL
    
    
    df=NULL
    df <- read.fst(paste0(dir_results,"results_Dispersal_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".fst"))

        for (pop in 1:npop){ # loop over popualtions
      
      demo <- NULL
      demo <- subset(df, Pop==pop)

      nyears <- max(demo$year,na.rm=TRUE)-1
      

      
      #-------------------------------#
      #  Composition des populations  #
      #-------------------------------#  
      
      
      for (i in 1:nyears){
        
        ## PARR 0+
 #       tmp<-subset(demo,Parr==1 & date==273 & AgeRiver<1 & year==i)
 #       Nparr0[i,pop] <- nrow(tmp)

        ## SMOLT 0+
 #       tmp<-subset(demo,Smolt==1 & date==90 & AgeRiver==1 & year==i)
 #       Nsmolt0[i,pop] <- nrow(tmp)

        ## ANADROMOUS
        tmp<-subset(demo,Returns==1 & date==273 & year==i)
        NRet[i,pop] <- nrow(tmp)
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
        
      } # end loop years
      
      tmp0[[pop]] <- Im.table
      
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
    
    
    # tmp1[[repmod]] <- NParr.table
    # tmp2[[repmod]] <- NSmolt.table
    # tmp3[[repmod]] <- list(Nhom=Nhomers.table, NIm=NIm.table, NEm=Reduce('+', tmp2), Im=tmp2)
    # tmp4[[repmod]] <- NRet.table

  } # end loop simul
  
  
#  nParr[[paste0(iEXPE)]] <- tmp1
#  nSmolt[[paste0(iEXPE)]] <- tmp2
  Mig[[paste0(iEXPE)]] <- tmp3
  nReturns[[paste0(iEXPE)]] <- tmp4
  
  
  #} # end scenario

### Save results
save(nReturns,Mig,file=paste0(dir_results,"DEMOGRAPHY",iEXPE,"_50.RData"))

q('no')
