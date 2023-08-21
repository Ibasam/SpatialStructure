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


#_________________ PARAMETERS _________________#

nSIMUL <- 50
nYear <- 40
nInit <- 10 # Nb years at initialization
npop <- 15


res_adult <- list() 

# SCENARIOS
iEXPE <- args

#for (iEXPE in EXPE){ # Loop over scenario 

  #_________________ DATA _________________#

  cat("Composition for ", iEXPE, "\n")
  
  for (repmod in 1:nSIMUL){ # loop over simulations
    cat("Simulation :",repmod," / ")
    perc <- (repmod/nSIMUL)*100
    if(perc %in% (seq(0,90,10))) cat(perc,"% ~ ","\n")
    if(perc == 100) cat(perc,"%","\n")
    
    df=NULL
    df <- read.fst(paste0(dir_results,"results_Dispersal_",iEXPE,"/Metapop_Sc",iEXPE,"_Sim",repmod,".fst"))
    
    res_adult_pop <- list()

    for (pop in 1:npop){
    
      # Remove initial values
      demo <- NULL
      demo <- subset(df, Pop==pop)

      id <- which(demo$Returns==1 & demo$date==273 & demo$Atsea==0) #1SW & MSW
      
      adult <- data.frame(ID = demo$ID[id]
                        , years = demo$year[id]
                        , Female = demo$Female[id]
                        , CollecID = demo$CollecID[id]
                        , gG = demo$gG[id]
                        , pG = demo$pG[id]
                        , motherID = demo$motherID[id]
                        , fatherID = demo$fatherID[id]
      )
      # MSW <- data.frame(ID = demo$ID[id2]
      #                   , years = demo$year[id2]
      #                   , Female = demo$Female[id2]
      #                   , Lf = demo$Lf[id2]
      #                   , gFmid3 = demo$gFmid3[id2]
      #                   , gFmid4 = demo$gFmid4[id2]
      #                   , CollecID = demo$CollecID[id2]
      #                   , AgeSea = demo$AgeSea[id2]
      #                   , gFmid1 = demo$gFmid1[id2]
      #                   , gFmid2 = demo$gFmid2[id2]
      #                   , gNeutral = demo$gNeutral[id2]
      #                   , gG = demo$gG[id2]
      #                   #, gG_sea = demo$gG_sea[id2]
      #                   #, gSLmid = demo$gSLmid[id2]
      #                   #, motherID = demo$motherID[id2]
      #                   #, fatherID = demo$fatherID[id2]
      #                   , pG = demo$pG[id2]
      #                   #, pG_sea = demo$pG_sea[id2]
      #                   #, W = demo$W[id2]
      #                   
      # )
      # 
      # SMOLT <- data.frame(ID = demo$ID[id3]
      #                     , years = demo$year[id3]
      #                     #, Female = demo$Female[id3]
      #                     , Lf = demo$Lf[id3]
      #                     #, AgeRiver = demo$AgeRiver[id3]
      #                     #, gG = demo$gG[id3]
      #                     #, gFmid1 = demo$gFmid1[id3]
      #                     #, gSLmid = demo$gSLmid[id3]
      #                     #, motherID = demo$motherID[id3]
      #                     #, fatherID = demo$fatherID[id3]
      #                     
      # )
      
      
      res_adult_pop[[pop]] <- adult
 #     res_MSW_pop[[pop]] <- MSW
 #     res_smolt_pop[[pop]] <- SMOLT
      
    } # end loop pop
    
    res_adult[[repmod]] <- res_adult_pop
 #   res_MSW[[repmod]] <- res_MSW_pop
 #   res_smolt[[repmod]] <- res_smolt_pop
    

  } # end loop simul
  

### Save results
save(res_adult,file=paste0(dir_results,"FITNESS",iEXPE,"_50.RData"))

q('no')
