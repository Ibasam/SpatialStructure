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


res_1SW <- list() 
res_MSW <- list() 
res_smolt <- list()

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
    
    res_1SW_pop <- list()
    res_MSW_pop <- list()
    res_smolt_pop <- list()

    for (pop in 1:npop){
    
      # Remove initial values
      demo <- NULL
      demo <- subset(df, Pop==pop)

      id <- which(demo$Returns==1 & demo$AgeSea<2 & demo$date==273 & demo$Atsea==0) #1SW
     # id2 <- which(demo$Returns==1 & demo$AgeSea>=2 & demo$date==273 & demo$Atsea==0) #MSW
 #     id3 <- which(demo$Smolt==1 & demo$date == 90) #smolt
      
      OSW <- data.frame(ID = demo$ID[id]
                        , years = demo$year[id]
                        , Female = demo$Female[id]
                        , Lf = demo$Lf[id]
                        , gFmid3 = demo$gFmid3[id]
                        , gFmid4 = demo$gFmid4[id]
                        , CollecID = demo$CollecID[id]
                        , AgeSea = demo$AgeSea[id]
                        , gFmid1 = demo$gFmid1[id]
                        , gFmid2 = demo$gFmid2[id]
                        , gNeutral = demo$gNeutral[id]
                        , gG = demo$gG[id]
                        #, gG_sea = demo$gG_sea[id]
                        #, gSLmid = demo$gSLmid[id]
                        #, motherID = demo$motherID[id]
                        #, fatherID = demo$fatherID[id]
                        , pG = demo$pG[id]
                        #, pG_sea = demo$pG_sea[id]
                        #, W = demo$W[id]
                        #, Fat = demo$Fat[id]

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
      #                   , W = demo$W[id2]
      #                   , Fat = demo$Fat[id2]
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
      
      
      res_1SW_pop[[pop]] <- OSW
 #     res_MSW_pop[[pop]] <- MSW
 #     res_smolt_pop[[pop]] <- SMOLT
      
    } # end loop pop
    
    res_1SW[[repmod]] <- res_1SW_pop
 #   res_MSW[[repmod]] <- res_MSW_pop
 #   res_smolt[[repmod]] <- res_smolt_pop
    

  } # end loop simul
  

### Save results
save(res_1SW,file=paste0(dir_results,"PHENOGENOTYPE",iEXPE,"_50.RData"))
#save(res_MSW,file=paste0(dir_results,"PHENOGENOTYPE",iEXPE,"_MSW_50.RData"))

q('no')
