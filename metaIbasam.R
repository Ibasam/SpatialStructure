#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
print(args)
if (length(args)==0){
  args <- ""
}

#_________________ PACKAGES _________________#
#install.packages("/metaIbasam_0.0.6.tar.gz", repos = NULL, type = "source")
library(metaIbasam)
library(doParallel)
library(MASS)
library(fst)
#_________________ LOAD IBASAM _________________#
source("Ibasam.R")

#_________________ SET RESULTS DIRECTORY _________________#
dir_results <- "results/"

#_________________ DATA & PARAMETERS _________________#
#RUN N SIMULATIONS
#nSIMUL=30
first = as.numeric(args[9])
last = as.numeric(args[10])
source("parIbasam.R")
#print(fishing.rates)

if(scenarioFishingRate>1 & sum(fishing.rates)==0){
  print("WARNINGS: Fishing rates should be > 0 !!!")
}


#Sc_name <-""
#Sc_name <- readline(prompt="Enter name: ")
# Sc_name <- paste0("Dispersal_",scenarioConnect,"_Fishing_",scenarioFishing,scenarioFishingRate,
#                   "_DiversityG_",scenarioDiversityG,"_DiversityTO_",scenarioDiversityTO,"_DivStruct_",
#                   scenarioDiversityStruct,"_Dist_",scenarioDist,"_Size_",scenarioSize)
Sc_name <- paste0("Dispersal_",scenarioConnect,scenarioFishing,scenarioFishingRate,
                  scenarioDiversityG,scenarioDiversityTO,
                  scenarioDiversityStruct,scenarioDist,scenarioSize)
if (Sc_name == "") {
  print('WARNING: no scenario name provided, temporary files and results will store into /tmp_ & /results_ folders')
  
}

# Create results folder
ifelse(!dir.exists(paste0(dir_results,'results_',Sc_name)),dir.create(paste0(dir_results,'results_',Sc_name)), FALSE)




### RUN


#for (sim in 1:nSIMUL) {
for (sim in first:last) {
    
  # Create temporary folder
  ifelse(dir.exists(paste0('tmp_',Sc_name)),unlink(paste0('tmp_',Sc_name), recursive = TRUE), FALSE)
  ifelse(!dir.exists(paste0('tmp_',Sc_name)),dir.create(paste0('tmp_',Sc_name)), FALSE)
  #system('mkdir tmp')
  
  #_________________ ENVIRONMENTAL CONDITIONS _________________#
  #source("code/river_climate_multi_Tair.R") #new file with the choice btw simulated data or observed data (external file) - but both based on air temperature
  if (Diversity_env==FALSE) {
    env <- env_new_false[[sim]]
  } else {
    env <- env_new_true[[sim]]
  }
  #source("code/river_climate_multi.R")
  
  #_________________ SIMULATIONS _________________#
  cl <- makeCluster(npop)#, outfile="log.txt")
  registerDoParallel(cl)
  
  results <- foreach(i=1:npop, .packages='metaIbasam', .verbose=T) %dopar% {
    
    tryCatch({
      
      Ibasam(nInit=nInit,nYears=nYears # Nb years simulated
             , npop = npop # Number of popualtions
             , Pop.o = i # Population of origin
             , rPROP = rPROP # Proportion de la taille de pop
             , CC.Temp=tempCC # Water Temperature increase (keep constant if 0)
             , CC.Amp=ampCC # Flow amplitude increase (keep constant if 1)
             , CC.Sea=seaCC #.75 # decreasing growth condition at sea (seqeunce from 1 to .75; keep constant if 1)
             , fisheries=fish.state, stage = fish.stage, fishing_rate=fishing.rates[i,] # fishing rates fn of population status
             , returning=TRUE # if TRUE, return R object contaning all data (very long)
             , plotting=FALSE,window=FALSE,success=FALSE,empty=TRUE
             , area=Area[i]
             , rmax= 10 #Piou #Rmax[14] #*M[i] #Rmax[i] #Rmax[15]*M[i]#10*M[i]     #Rmax reference Scorff * coeff multiplicateur
             , alpha= 0.1 #Piou #*M[i] #Alpha[i] #changé par 0.1 pour laisser comme papier Piou 2012 #alpha[15]*M[i]   #alpha reference Scorff * coeff multiplicateur
             , pstray=pstray[i,]  # Dispersal probability #ligne popo
             , propG = propG[i]   #Genetic change at initialization - al - 12/02/2021
             , propTO = propTO[i]   #Trade-off parameter change at initialization - al - 08/12/2021
      )
      
    }, error = function(e) return(paste0("The population '", i, "'",  " caused the error: '", e, "'")))
    #cat(dput(i), file = paste0("debug_file_", i, ".txt"))
  }
  #save(results,file=paste0(dir_results,"results_",Sc_name,"/Metapop_Sc",scenarioConnect,scenarioFishing,scenarioFishingRate,"_Sim",sim,".RData"))
  #save(results,file=paste0("results/Metapop_Sc",scenarioConnect,scenarioEnvi,scenarioFishing,"rhoT",rhoT,"rhoF",rhoF,"_Sim",10,"_observed.RData"))
  
  # Extract pop results from list
  res.pop=NULL
  #fishing_rates=NULL
  for (p in 1:npop){
   res.pop[[p]] <- results[[p]]$pop
   res.pop[[p]]$Pop <- p

   # fishing_tmp<-NULL
   # fishing_tmp <- results[[p]]$fishing_rate
   # fishing_rates <- rbind(fishing_rates, fishing_tmp)
  }
  # binding columns together
  df <- do.call(rbind, res.pop)
  # converting to a dataframe
  df <- as.data.frame(df)
  # Save to fst format 
  compress_rate=100 #0: no compress / 100: max
  write.fst(df,paste0(dir_results,"results_",Sc_name,"/Metapop_Sc",scenarioConnect,scenarioFishing,scenarioFishingRate,scenarioDiversityG,scenarioDiversityTO,scenarioDiversityStruct,scenarioDist,scenarioSize,"_Sim",sim,".fst"),compress_rate)
  #save(fishing.rates,file=paste0(dir_results,"results_",Sc_name,"/Fishing_rates_Sc",scenarioConnect,scenarioFishing,scenarioFishingRate,scenarioDiversityG,scenarioDiversityTO,"_Sim",sim,".RData"))
  

  
  # Cleaning
  stopCluster(cl)
  rm(cl)
  gc() # clean memory
  system(paste0('rm -R tmp_',Sc_name)) #if want to not remove the folder tmp
  
}

#if (args > 1) { q('no') } # close R session 
q('no')
