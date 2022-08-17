Ibasam <-
  function (nInit #nbr d'année pour faire burn-in
            ,nYears #nbr d'année pour faire simulation
            , npop # Number of populations
            , Pop.o # Population of origin
            , rPROP
            #, dataEnv = TRUE
            , CC.Temp # effects of climate change in freshwater on temperature
            , CC.Amp # effects of climate change in freshwater on flow
            , CC.Sea # effects of climate change at sea
            , fisheries = TRUE, stage = TRUE, fishing_rate = fishing.rates ## fishing effects
            , plotting = TRUE, window = FALSE, returning = TRUE, success = FALSE, empty = TRUE
            , area 
            , rmax
            , alpha
            , pstray # Dispersal probability
            , propG #genetic change at initialization # al - 12/02/2021
            , propTO #trade-off parameter change at initialization # al - 08/12/2021
  ) 
  {
    
    #Initialization & Preparation:
    empty()
    def <- defaultParameters()
    
    #### PARAMETERS ####
    # If maxRIV <0 or maxSEA <0, NO trade-offs
    kappaRIV=0.001
    kappaSEA=0.001
    
    heriRIV=0.14
    heriSEA=0.14
    
    maxSEA=50
    sigSEA=100
    
    maxRIV=5 #-1
    sigRIV=3.7 * (1+ propTO) #edit al TO - 09-03-21
    
    SP0=0.9841606*1.0025 # edit buoro -> suvival too high and density of 0+ parr too high with 1,005
    SP1=0.9914398*1.0025 # edit buoro
    SP1S=0.9967923*1.002
    SP1M=0.9863295*1.002
    SPnM=0.9911798*1.002
    SPn=0.99775*1.002
    
    #les paramètres qu'il faut changer 
    #aire de production * proportion de l'aire qu'on prend rPROP
    def$envParam[9] <- area*rPROP
    #survies
    def$colParam[13:18] <- c(SP0,SP1,SP1S,SP1M,SPnM,SPn) #daily survival prob parr 0 to 0.5; daily survival prob parr 0.5 to 1.0; daily survival prob smolts 0.5 to run; daily survival prob parr mature 0.5 to 1.0; daily survival prob parr mature n.5 to n+1; daily survival prob parr any other situations
    #19:24 = forme du compromis (maxRIV,sigRIV,kappaRIV,maxSEA,sigSEA,kappaSEA)
    def$colParam[19:24] <- c(maxRIV,sigRIV,kappaRIV,maxSEA,sigSEA,kappaSEA)
    #alpha
    def$colParam[40] <- alpha
    #Rmax: Neggsmax
    def$colParam[41] <- area*rPROP*rmax #def$envParam[9]*2   #2 oeufs / m2
    #63 = h?ritabilit? de la croissance en rivi?re
    def$colParam[63] <- heriRIV
    #67 = h?ritabilit? de la croissance en mer
    def$colParam[67] <- heriSEA
    
    
    # Genetic traits at initialization - al - 12/02/2021
    def$colParam[61] <- log(exp(def$colParam[61]) * (1+ propG))#gMean # Growth in river ; gMean
    def$colParam[65] <- log(exp(def$colParam[65]) * (1+ propG))#gMean # Growth at sea ; gGmeanSea
    def$colParam[69] <- def$colParam[69] * (1+ propG)#LmidSm           ->Norm reaction smoltification
    def$colParam[77] <- def$colParam[77] * (1+ propG)#parrMalesTreshold #Fmid[0] threshold of parr males
    def$colParam[80] <- def$colParam[80] * (1+ propG)#parrFemalesTreshold #Fmid[1] threshold of parr females
    def$colParam[83] <- def$colParam[83] * (1+ propG)#anadMalesTreshold #Fmid[2] threshold of salmon males
    def$colParam[86] <- def$colParam[86] * (1+ propG)#anadFemalesTreshold #Fmid[3] threshold of salmon females
    
    
    #on retrouve l'aire et la proportion d'aire un peu partout, alors il faut aussi les changer
    #c'est en fait ce qui a été fait ici
    def$gParam[1] <- round(area*rPROP*8*0.15)
    def$parrParam[1] <- round(area*rPROP*8*0.15*0.011)
    def$smoltsParam[1] <- round(area*rPROP*8*0.15*0.03)
    def$grilseParam[1] <- round(area*rPROP*8*0.15*0.003)
    def$mswParam[1] <- round(area*rPROP*8*0.15*0.0005)
    
    #### ENVIRONMENT ####
    #mm <- river_climate_model(nInit + nYears + 1, CC.Temp, CC.Amp) #script with no choice of environmental covariation (not forget to commence in metaIbasam.R the source "river_climate_multi")
    # if (dataEnv==TRUE){
    #   mm <- river_climate_model(nInit + nYears + 1, CC.Temp, CC.Amp)
    # } else {
    #   mm <- river_climate_model(nInit + nYears + 1, CC.Temp, CC.Amp)
    # }
    # load("data/river_climate_multi.RData") #edit alamarins for the script with choice of env covariation BUT NO NEEDED (source in metaIbasam.R)
     mm=NULL #edit alamarins for the script with choice of env covariation
     mm$flow <- env$flow[,Pop.o] #edit alamarins for the script with choice of env covariation #now in natural scale (metaibasam v6)
     mm$temperatures <- env$temperatures[,Pop.o] #edit alamarins for the script with choice of env covariation
    
     mm$module <- env$module[Pop.o]
     def$envParam[7] <- 0.2 * mm$module # Critical_RelativeFlow
     def$envParam[15] <- 0.1 * mm$module #0 * mm$module # CritInfFlow_
     def$envParam[16] <- 7 * mm$module #20 * mm$module #CritSupFlow_
     
     #aRt
     #def$envParam[10] <- 0.0001936 #al - 06/10/2020
     
     #dr
     def$envParam[3] <- 0.5 #* (1- prop) #edit al - 18/03/21 opposite for local adap #0.427 #al - 06/10/2020
     
     
     # Oceanic growth conditions:
     MeanNoiseSea <- c(rep(1,nInit),seq(1,CC.Sea,length=nYears))
     mm$MeanNoiseSea <- MeanNoiseSea
     
     # 
    Reset_environment()
    Prepare_environment_vectors(mm$temperatures, mm$flow)
    setup_environment_parameters(def$envParam)
    setup_collection_parameters(def$colParam)
    
 
    
    ## Define fishing rates  
    if (fisheries) {
      if(stage){
        rates <- cbind(
          grilses=c(rep(fishing_rate[1],nInit),rep(fishing_rate[2],nYears))
          ,msw=c(rep(fishing_rate[3],nInit),rep(fishing_rate[4],nYears)) #(tbeauchard 13/04/2021 : ajout du taux d'exploitation initial et au bout de 10ans)
        )         
      } else {
        rates <- cbind(
          Small=c(rep(fishing_rate[1],nInit),rep(fishing_rate[1],nYears))
          ,Med=c(rep(fishing_rate[2],nInit),rep(fishing_rate[2],nYears))
          ,Big=c(rep(fishing_rate[3],nInit),rep(fishing_rate[3],nYears))
        )
      }
    }
    
    #### INITIALIZING POPULATION ####
    set_collecID(Pop.o) # provide ID number for each popualtion (variable CollecID) to avoid duplicate individulas ID with migrants
    time_tick(90)
    add_individuals(def$gParam)
    add_individuals(def$parrParam)
    add_individuals(def$smoltsParam)
    go_summer()
    popo <- observe()
    add_individuals(def$grilseParam)
    add_individuals(def$mswParam)
    go_winter()
    
    popa <- observe()
    popb <- observe_redds()# edit mbuoro
    #source("code/memory_success.R") #edit alamarins
    #popc <- memory_success()#edit alamarins
    if (returning || success) {
      results <- popa #before: observe()
      redds <- popb# edit mbuoro #before:observe_redds()
    #  emerg <- popc#edit alamarins
    }
    
    
    ratios <- matrix(NA, nrow = nInit+nYears, ncol = 4)
    winterM <- matrix(NA, nrow = nInit+nYears, ncol = 6)
    summerM <- matrix(NA, nrow = nInit+nYears, ncol = 18)
    ally <- summarize.oneyear(popo, popa)
    sptm <- NULL
    
    
    ## RUN
    pb   <- txtProgressBar(1, nYears+nInit, style=3) # initilazing progress bar
    N<-NULL #
    for (y in 1:(nYears+nInit)) {
      #cat("Year: ",y,"of ",nYears, "\n")
      setTxtProgressBar(pb, y) # progress bar
      
      # Oceanic growth conditions:
      def$envParam[1] <- MeanNoiseSea[y]
      setup_environment_parameters(def$envParam)
      
      ptm <- proc.time()
      spring()
      #emerg <- observe() #edit alamarins to get nbr of emergents (pb SRR to solve)
      summer()
      
      #### FISHING ####
      popo <- observe() # state BEFORE fisheries
      if (fisheries) {
        fishing(rates[y,])
      }
      
      #popo <- observe() # state AFTER fisheries        
      if (returning || success) {
        results <- rbind(results, popo)
        #results <- rbind(results, emerg) #edit alamarins to get nbr of emergents (pb SRR to solve)
        
      }
      
      ratios[y, ] <- unlist(proportions.population(popo))
      summerM[y, ] <- unlist(important.indicator.summer.population(popo))
      
      autumn()
      winter()
      
      #trying to understand the model bug (stop run when few nbr of individuals which stops the run of all pops)
      beug<- observe() #edit alamarins
      n<-nrow(beug) #edit alamarins 
      N<-c(N,n) #edit alamarins
      tmp<-cbind(1:y,N) #edit alamarins
      write.table(tmp, file = paste0("tmp_",Sc_name,"/BeugPop_",Pop.o,".txt") , sep = "\t", 
                  row.names = TRUE, col.names = NA) #edit alamarins to get the temporal evolution of nbr of ind
      # 
      # 
      #### STRAYING ####
      #emmigrants("nom de fichier", straying_rates for 1SW & MSW)
      #pause("nom de fichier")
      #immigrants("nom de fichier")
      # emmigrants(paste0("tmp/Pop_",Pop,"_",y,".txt"),c(0.1,0.1))
      # popb <- observe()
      # pause(paste0("tmp/Pop_",Pop_im,"_",y,".txt")) # R script to pause the execution of Ibasam until immigrant file (e.g. mig_AtoB) is created in a specific folder
      # immigrants(paste0("tmp/Pop_",Pop_im,"_",y,".txt"))
      # popc <- observe()
      
      for (Pop.e in 1:npop){
        # Pop.o: population of origin
        # Pop.e: emigrate to population Pop.e
        if(Pop.e == Pop.o) { 
          next 
        } else {
          ## boucle if edit alamarins to force the model to create the file even if few ind in a pop and bug pop, the others pop will continue the simul
          # if(n<=5) {
          #   #file.create(file = paste0("tmp/Mig_",Pop.o,"-",Pop.e,"_",y:(nYears+nInit),".txt"))            #i c'est l'indice pour nSimu définit dans scriptR
          #   emfile <- paste0("tmp/Mig_",Pop.o,"-",Pop.e,"_",y,".txt") 
          #   
          #   file.create(file = paste0("tmp/Crash_",Pop.o,"-",Pop.e,"_",y,".txt"))            #i c'est l'indice pour nSimu définit dans scriptR 
          #   write.table(beug, file = paste0("tmp/Crashwinter_",Pop.o,".csv") , sep = "\t",
          #               row.names = TRUE, col.names = NA)
          #   write.table(popo, file = paste0("tmp/Crashsummer_",Pop.o,".csv") , sep = "\t",
          #               row.names = TRUE, col.names = NA)
          # } else {
          emfile <- paste0("tmp_",Sc_name,"/Mig_",Pop.o,"-",Pop.e,"_",y,".txt") 
          #i c'est l'indice pour nSimu définit dans scriptR 
          #y c'est le numéro de l'année en incluant celles de burnin 
          #pstray <- c(0.1,0.1)
          if (y<= nInit) {
            emmigrants(emfile,0) #no dispersal during initialization phase #edit mbuoro
          } else {
            emmigrants(emfile, pstray[Pop.e])
          }
          #}
          
        } # end if
      } # end Pop.e
      
      #pope <- observe()
      
      for (Pop.i in 1:npop){
        # Pop.o: population of origin
        # Pop.i: immigrate from population Pop.i
        if(Pop.i == Pop.o) { 
          next 
        } else {
          imfile <- paste0("tmp_",Sc_name,"/Mig_",Pop.i,"-",Pop.o,"_",y,".txt")
          pause(imfile) # R script to pause the execution of Ibasam until immigrant file (e.g. mig_AtoB) is created in a specific folder
          immigrants(imfile)
        } # end if
      } # end Pop.i
      
      # popi <- observe()
      # if (returning || success) {
      #   results <- rbind(results, pope)
      # }
      
      
      popa <- observe() 
      popb <- observe_redds()# edit mbuoro
     # popc <- memory_success()#edit alamarins
      if (returning || success) {
        results <- rbind(results, popa)
        redds <- rbind(redds, popb)# edit mbuoro
       # emerg <- rbind(emerg, popc)#edit alamarins
      }
      
      winterM[y, ] <- unlist(important.indicator.winter.population(popa))
      ally <- append.oneyear(popo, popa, ally)
      sptm <- rbind(sptm, proc.time() - ptm)
    }
    
    #### PLOT ####
    if (plotting) {
      pdf(paste('tmp/Res_Pop',Pop.o,'.pdf',sep=''))
      # ecrase à chaque fois l'ancier pdf, ne garde que celui de la dernière simulation
      op <- par(mfrow = c(2, 2))
      plot_proportions_population(ratios, window = window)
      plot_winterM(winterM, window = window)
      plot_summerM(summerM, window = window)
      plotevolution(ally, window = window)
      par(mfrow = c(2, 1))
      if (success) {
        newwindow(window)
        suc <- temporal_analyse_origins(results, 1:nYears, 
                                        plotting = plotting, titles = "Strategy success through time")
      }
      newwindow(window)
      plot(ts(sptm[, 1]), main = "CPU time needed per year", 
           ylab = "seconds", xlab = "years", bty = "l", sub = paste("Total:", 
                                                                    round(sum(sptm[, 1]), 3)))
      lines(lowess(sptm[, 1]), col = 2, lty = 2)
      par(op)
      dev.off()
    }
    if (returning) {
      #return(list(results,mm))
      #return(list("pop"=results,"redds"=redds,"env"=mm,"alpha"=def$colParam[40],"Rmax"=def$colParam[41])) # edit mbuoro #edit alamarins emerg
      return(list("pop"=results))
    } else {
      invisible(NULL)
    }
  }
