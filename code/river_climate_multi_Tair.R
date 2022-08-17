#param relation Tair-Twater_flow Scorff (adjusted on data from Nicolas Jeannot and Guillaume Thirel data on Scorff river)
#-> some years with missing data -> data reconstructed through hierarchical bayesian model)
#AR order 15
#Q is in natural scale
#load(paste0("data/param_logistiq_relation_airtemp_waterflow_watertemp.Rdata"))

## median of parameters over 10000 iterations of one markov chain ##
theta1 <- 0.8686904 #median(fit.mcmc[[1]][10000:20000,"theta1"])
theta2 <- 0.1782612 #median(fit.mcmc[[1]][10000:20000,"theta2"])
theta3 <- 2.376725 #median(fit.mcmc[[1]][10000:20000,"theta3"])
alpha1 <- 28.08205 #median(fit.mcmc[[1]][10000:20000,"alpha1"])
alpha2 <- 0.1867286 #median(fit.mcmc[[1]][10000:20000,"alpha2"])
#alpha3 <- fit.mcmc[[1]][smp,"alpha3"]
sd.theta <- 0.7132355 #median(fit.mcmc[[1]][10000:20000,"sd.theta"])
sd.alpha <- 1.217471 #median(fit.mcmc[[1]][10000:20000,"sd.alpha"])

# a<-NULL
# for (i in 1:16) {
#   a[i] <- median(fit.mcmc[[1]][10000:20000,paste0("a[",i,"]")])
# }
a <- c(0.03336559, 0.05000371, 0.02577071, 0.01599289, 0.01060309, 0.005616591, 0.005241847, 0.004025856, 0.003048239, 0.002696868, 
     0.003689901, 0.002452315, 0.001433584, 0.002061772, 0.0009726248, 0.005033227)
b <- -0.2315342 #median(fit.mcmc[[1]][10000:20000,"b"])
c <- -2.153799 #median(fit.mcmc[[1]][10000:20000,"c"])
#d <- fit.mcmc[[1]][smp,"d"]
sigmaTwater <- 0.8409419 #median(fit.mcmc[[1]][10000:20000,"sd"])



if (Obs==FALSE) {
  
  #_________________ FUNCTION _________________#
  sinus_model_resid_ar_multi <-  function (m, amp, phase, ar, eps, day, nday) 
  {
    error <- numeric(length(day))
    tm <- m + amp * sin(2 * pi * (day - phase)/nday)
    error[1] <- 0
    for (dd in day[-1]) {
      error[dd] <- ar[1] * error[dd - 1] + eps[dd]
      tm[dd] <- tm[dd] + error[dd]
    }
    return(tm)
  }
  
  
  ### PARAMETERS####
  #Rivers parameters Tair (adjusted on data from G Thirel)
  pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Ellee", "Scorff", "Blavet")
  
  #mTair <- ampTair <- phaseTair <- arTair <- sigmaTair <- NULL
  if (Diversity_env == TRUE) {
    #for (pop in 1:15) {
      #load(paste0("data/param_airtemperature_single_constant_AR_nimble_",pops[pop],".Rdata"))
      
      mTair <- c(10.85441, 10.46192, 10.77484, 10.20836, 10.61513, 10.40890, 10.57859, 10.75085, 10.30889, 11.43305, 10.83083, 10.96849, 10.58954, 10.70966, 10.29672)#c(mTair,median(samplesList[[1]][5000:10000,"alpha0"]))
      ampTair <- c(5.625823, 5.681926, 5.467652, 5.562698, 5.405622, 5.190783, 5.175363, 5.161347, 5.584885, 5.254536, 5.477578, 5.481080, 5.788681, 5.869587, 5.792893) #c(ampTair,median(samplesList[[1]][5000:10000,"amp0"]))
      phaseTair <- c(120.1278, 119.3717, 120.7393, 119.6634, 120.6235, 121.0866, 121.2506, 121.3283, 119.1253, 121.0813, 119.2211, 119.3495, 117.9483, 117.5118, 118.6635) #c(phaseTair,median(samplesList[[1]][5000:10000,"phase0"]))
      arTair <- c(0.7871263, 0.7924448, 0.7831106, 0.7948227, 0.7871766, 0.7857606, 0.7861066, 0.7846760, 0.7949550, 0.7997064, 0.7925223, 0.7948378, 0.7925832, 0.7895524, 0.7944258) #c(arTair,median(samplesList[[1]][5000:10000,"ar1"]))
      sigmaTair <- c(1.710117, 1.706664, 1.690047, 1.656230, 1.655227, 1.605951, 1.595978, 1.599738, 1.685187, 1.517897, 1.644425, 1.626768, 1.731746, 1.779085, 1.727437) #c(sigmaTair,median(samplesList[[1]][5000:10000,"sd"]))
    #}
    # for (pop in c(2,4, 8, 12)) { ### test better conditions on some rivers #al-10-11-20
    #   mTair[pop] <- mTair[pop] + 2
    # }
  } else {
    #load(paste0("data/param_airtemperature_single_constant_AR_nimble_",pops[14],".Rdata"))
    
    mTair <- rep(10.70966,npop) #median(samplesList[[1]][5000:10000,"alpha0"])
    ampTair <- rep(5.869587,npop) #median(samplesList[[1]][5000:10000,"amp0"])
    phaseTair <- rep(117.5118,npop) #median(samplesList[[1]][5000:10000,"phase0"])
    arTair <- rep(0.7895524,npop) #median(samplesList[[1]][5000:10000,"ar1"])
    sigmaTair <- rep(1.779085,npop) #median(samplesList[[1]][5000:10000,"sd"])
  }

  
  #Rivers parameters Flow:(adjusted on data from G. Thirel) -> log(Q/module)
  pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Ellee", "Scorff", "Blavet")
  
  #mF <- ampF <- phaseF <- arF <- sigmaF <- NULL
  if (Diversity_env == TRUE) {
    #for (pop in 1:15) {
      #load(paste0("data/param_flow_single_constant_AR_nimble_",pops[pop],".Rdata"))
      mF <- c(-0.6011924, -0.5180006, -0.6073302, -0.4902582, -0.3775791, -0.3155432, -0.4559050, -0.3763208, -0.8321171, -0.4867163, -0.5639587, -0.3998210, -0.5010848, -0.4241705, -0.5170546) #c(mF, median(samplesList[[1]][5000:10000,"alpha0"]))
      ampF <- c(1.1877145, 1.1359234, 1.1795817, 1.1145294, 0.9669853, 0.7620662, 1.0666736, 0.8742592, 1.5267847, 1.1233765, 1.2190544, 0.9777526, 1.1045493, 0.9943968, 1.1791851) #c(ampF, median(samplesList[[1]][5000:10000,"amp0"]))
      phaseF <- c(329.4959, 327.6700, 325.3111, 321.9974, 326.0787, 323.6536, 320.4378, 312.7823, 311.7362, 321.0212, 318.4419, 317.4873, 323.0439, 325.9451, 323.0537) #c(phaseF, median(samplesList[[1]][5000:10000,"phase0"]))
      arF <- c(0.9496850, 0.9384839, 0.9221304, 0.9246938, 0.9499700, 0.8992473, 0.9464243, 0.9257478, 0.9328857, 0.9460142, 0.9477278, 0.9307698, 0.9532685, 0.9391304, 0.8986100) #c(arF, median(samplesList[[1]][5000:10000,"ar1"]))
      sigmaF <- c(0.2241899, 0.2236890, 0.2745719, 0.2404288, 0.1759966, 0.2269506, 0.2027411, 0.2237388, 0.3069457, 0.2040074, 0.2208548, 0.2113761, 0.1991157, 0.2100546, 0.2950159) #c(sigmaF, median(samplesList[[1]][5000:10000,"sd"]))
    #}
  } else {
    #load(paste0("data/param_flow_single_constant_AR_nimble_",pops[14],".Rdata"))
    mF <- rep(-0.4241705,npop) #median(samplesList[[1]][5000:10000,"alpha0"])
    ampF <- rep(0.9943968,npop) #median(samplesList[[1]][5000:10000,"amp0"])
    phaseF <- rep(325.9451,npop) #median(samplesList[[1]][5000:10000,"phase0"])
    arF <- rep(0.9391304,npop) #median(samplesList[[1]][5000:10000,"ar1"])
    sigmaF <- rep(0.2100546,npop) #median(samplesList[[1]][5000:10000,"sd"])
  }

  if (Diversity_env == TRUE) {
    module <- c(2.7154143,  5.1510163 , 1.7321199 , 4.7378904 , 0.8144686 , 0.4467169 , 2.8034946 , 5.9937238 , 2.2431318 , 1.4613574 , 4.9348077 , 3.6451468 , 9.6662743, 4.9512308, 10.1685114) #vector of module or rivers (from data G Thirel)
  } else {
    module <- c(rep(4.9512308,npop)) #Scorff
  }
  
  
  
  #choice of the covariance between rivers conditions (through the residuals epsilon)
  covMatTair <- sigmaTair %*% t(sigmaTair) * corMat
  #covMatTair
  covMatF <- sigmaF %*% t(sigmaF) * corMat
  #covMatF
  
  eTair <- mvrnorm(365*(nInit + nYears + 1),mu = rep(0,npop),Sigma = covMatTair, empirical = TRUE)
  eF <- mvrnorm(365*(nInit + nYears + 1),mu = rep(0,npop),Sigma = covMatF, empirical = TRUE)
  
  #eTair <- rnorm(365*(nInit + nYears + 1),0, sigmaTair)
  
  
  
  
  ##function to simulate environmental conditions (log(Q/module) & Tair)
  
  
  river_climate_model_multi <- function (npop, nInit, nYears, CC_Temp, CC_Amp) 
  {
    
    nyear=nInit + nYears + 1
    nday <- 365
    day <- 1:(365 * nyear) #rep(1:365,nyear)

    temperatures=NULL
    flow=NULL
    
    #compute Tair for each river (data G Thirel))
    for (pop in 1:npop){
      
    if(missing(CC_Temp)) {
      mTair <- mTair[pop]
    } else {
      #mT<-seq(mT[pop],mT[pop]+ CC_Temp,length=365*nyear)  #MeanT
      mTair<-c(
        seq(mTair[pop],mTair[pop],length=365*nInit)  #MeanT
        ,seq(mTair[pop],mTair[pop]+ CC_Temp,length=365*(nYears + 1))  #MeanT
      )
    }
    tair <- sinus_model_resid_ar_multi(mTair, ampTair[pop], phaseTair[pop], arTair[pop], eTair[,pop], day, nday)
    tair_year <- data.frame(cbind(tair, rep(1:nyear, each=365)))
    colnames(tair_year) <- c("tair","year")
    
    #compute logrelflow for each river (data G Thirel) 
      if(missing(CC_Amp)) {
        ampF <- ampF[pop]
      } else {
        #mF <- seq(mF[pop],1.010139* CC_Amp,length=365*nYears)
        #  ampF <- seq(ampF[pop],ampF[pop]* CC_Amp,length=365*nyear)
        ampF <- c(
          seq(ampF[pop],ampF[pop],length=365*nInit)
          ,seq(ampF[pop],ampF[pop]* CC_Amp,length=365*(nYears + 1))
        )
      }
      #csF <- 327.836285
      #alpF <- 0.964034622660953
      #betF <- 0.118343974744226
     
      
      tmplogrelflow <- sinus_model_resid_ar_multi(mF[pop], ampF, phaseF[pop], arF[pop], eF[,pop], day, nday)
      
      ##transform logrelflow in flow (natural scale)
      #if (Diversity_env == TRUE) {
        tmpflow <- exp(tmplogrelflow) * module[pop]
      #} else {
        #tmpflow <- exp(tmplogrelflow) * module
      #}
 
      flow_year <- data.frame(cbind(tmpflow, rep(1:nyear, each=365)))
      colnames(flow_year) <- c("flow","year")
      
      ##compute IQ for water temperature computing
      #if (Diversity_env == TRUE) {
        IQ <- 1/((module[pop]/10)+tmpflow)
      #} else {
        #IQ <- 1/((module/10)+tmpflow)
      #}
      
      
      
      ##compute water temperature

      min_Tair_year=(aggregate(tair~year,tair_year, min, na.rm=T)$tair)-mean(aggregate(tair~year,tair_year, min, na.rm=T)$tair)
      max_Tair_year=(aggregate(tair~year,tair_year, max, na.rm=T)$tair)-mean(aggregate(tair~year,tair_year, max, na.rm=T)$tair)
      min_flow_year=(aggregate(flow~year,flow_year, min, na.rm=T)$flow)-mean(aggregate(flow~year,flow_year, min, na.rm=T)$flow)
      #max_flow_year=(aggregate(flow~year,rivers[which(rivers$river==pop),], max, na.rm=T)$flow)-mean(aggregate(flow~year,rivers[which(rivers$river==pop),], max, na.rm=T)$flow)
          
      min_Tair=max_Tair=min_flow=max_flow=NULL
      for (i in 1:length(min_Tair_year)) {
          min_Tair=c(min_Tair, rep(min_Tair_year[i],365))
          max_Tair=c(max_Tair, rep(max_Tair_year[i],365))
          min_flow=c(min_flow, rep(min_flow_year[i],365))
          #max_flow=c(max_flow, rep(max_flow_year[i],365))
      }
          
      #min
      mu.theta <- theta1 + (theta2 * min_Tair) + (theta3 * min_flow)
      theta <- mu.theta + rnorm(1,0,sd.theta) #big source of variability btw simulations ?
      #max
      mu.alpha <- alpha1 + (alpha2 * max_Tair) #+ (alpha3 * min_flow)
      alpha <- mu.alpha + rnorm(1,0,sd.alpha)
          
         
      #compute water temperature
      x=NULL
      x[1] <- sum(a[1] * tair[1]) + b * IQ[1] + c
      x[2] <- sum(c(a[1] * tair[2], a[2] * tair[1])) + b * IQ[2] + c
      x[3] <- sum(c(a[1] * tair[3], a[2] * tair[2], a[3] * tair[1])) + b * IQ[3] + c
      x[4] <- sum(c(a[1] * tair[4], a[2] * tair[3], a[3] * tair[2], a[4] * tair[1])) + b * IQ[4] + c
      x[5] <- sum(c(a[1] * tair[5], a[2] * tair[4], a[3] * tair[3], a[4] * tair[2], a[5] * tair[1])) + b * IQ[5] + c
      x[6] <- sum(c(a[1] * tair[6], a[2] * tair[5], a[3] * tair[4], a[4] * tair[3], a[5] * tair[2], a[6] * tair[1])) + b * IQ[6] + c
      x[7] <- sum(c(a[1] * tair[7], a[2] * tair[6], a[3] * tair[5], a[4] * tair[4], a[5] * tair[3], a[6] * tair[2], a[7] * tair[1])) + b * IQ[7] + c
      x[8] <- sum(c(a[1] * tair[8], a[2] * tair[7], a[3] * tair[6], a[4] * tair[5], a[5] * tair[4], a[6] * tair[3], a[7] * tair[2], a[8] * tair[1])) + b * IQ[8] + c
      x[9] <- sum(c(a[1] * tair[9], a[2] * tair[8], a[3] * tair[7], a[4] * tair[6], a[5] * tair[5], a[6] * tair[4], a[7] * tair[3], a[8] * tair[2], a[9] * tair[1])) + b * IQ[9] + c
      x[10] <- sum(c(a[1] * tair[10], a[2] * tair[9], a[3] * tair[8], a[4] * tair[7], a[5] * tair[6], a[6] * tair[5], a[7] * tair[4], a[8] * tair[3], a[9] * tair[2], a[10] * tair[1])) + b * IQ[10] + c
      x[11] <- sum(c(a[1] * tair[11], a[2] * tair[10], a[3] * tair[9], a[4] * tair[8], a[5] * tair[7], a[6] * tair[6], a[7] * tair[5], a[8] * tair[4], a[9] * tair[3], a[10] * tair[2], a[11] * tair[1])) + b * IQ[11] + c
      x[12] <- sum(c(a[1] * tair[12], a[2] * tair[11], a[3] * tair[10], a[4] * tair[9], a[5] * tair[8], a[6] * tair[7], a[7] * tair[6], a[8] * tair[5], a[9] * tair[4], a[10] * tair[3], a[11] * tair[2], a[12] * tair[1])) + b * IQ[12] + c
      x[13] <- sum(c(a[1] * tair[13], a[2] * tair[12], a[3] * tair[11], a[4] * tair[10], a[5] * tair[9], a[6] * tair[8], a[7] * tair[7], a[8] * tair[6], a[9] * tair[5], a[10] * tair[4], a[11] * tair[3], a[12] * tair[2], a[13] * tair[1])) + b * IQ[13] + c
      x[14] <- sum(c(a[1] * tair[14], a[2] * tair[13], a[3] * tair[12], a[4] * tair[11], a[5] * tair[10], a[6] * tair[9], a[7] * tair[8], a[8] * tair[7], a[9] * tair[6], a[10] * tair[5], a[11] * tair[4], a[12] * tair[3], a[13] * tair[2], a[14] * tair[1])) + b * IQ[14] + c
      x[15] <- sum(c(a[1] * tair[15], a[2] * tair[14], a[3] * tair[13], a[4] * tair[12], a[5] * tair[11], a[6] * tair[10], a[7] * tair[9], a[8] * tair[8], a[9] * tair[7], a[10] * tair[6], a[11] * tair[5], a[12] * tair[4], a[13] * tair[3], a[14] * tair[2], a[15] * tair[1])) + b * IQ[15] + c
          
      muT=tmpT=NULL
      for (t in 1:15) {
        muT[t] <- theta[t] + ((alpha[t] - theta[t])/(1 + exp(-x[t])))
        tmpT[t] <- muT[t] + rnorm(1,0,sigmaTwater)
      }
          
          
      for (t in 16:length(day)) {
        x[t] <- sum(c(a[1] * tair[t], a[2] * tair[t-1], a[3] * tair[t-2], a[4] * tair[t-3], a[5] * tair[t-4], a[6] * tair[t-5], a[7] * tair[t-6], a[8] * tair[t-7], a[9] * tair[t-8], a[10] * tair[t-9], a[11] * tair[t-10], a[12] * tair[t-11], a[13] * tair[t-12], a[14] * tair[t-13], a[15] * tair[t-14], a[16] * tair[t-15])) + b * IQ[t] + c
        muT[t] <- theta[t] + ((alpha[t] - theta[t])/(1 + exp(-x[t])))
        tmpT[t] <- muT[t] + rnorm(1,0,sigmaTwater)
      }
    
      
      tmpT[tmpT < 0] <- 0.01
      tmpflow[tmpflow <0] <- 0
      # tmplogrelflow <- tmplogrelflow -  mF[pop] #correction due to values not centered on 0
      # tmplogrelflow[tmplogrelflow > 3] <- 3
      # tmplogrelflow[tmplogrelflow < -3] <- -3
      # 
      # if (pop %in% c(2,4, 8, 12)) { ### test better conditions on some rivers #al-10-11-20
      #   tmpT <- tmpT + 2
      # }
      
      temperatures <- cbind(temperatures,tmpT)    
      flow <- cbind(flow,tmpflow)
      
    } # end loop pop
    
    return(list(mTair=mTair,ampF=ampF,temperatures = temperatures, flow = flow, module=module))
    
} #end function river_climate_model_multi

  env <- river_climate_model_multi(npop, nInit, nYears, tempCC, ampCC) 
  
} else { #end if Obs=FALSE   

  
  ### DATA####
  ###choose the model GR
  ###load simulated data for each river (keep the same order of pop (Leff-> Blavet))
  pops <- c(
    "Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze",
    "Elorn", "Aulne", "Goyen", "Odet",
    "Aven","Ellee", "Scorff", "Blavet"
  )
  code <- c("J1813010","J1721720","J2023010","J2233010","J2314910","J2404010","J2723010","J3413030","J3601810","J4014010","J4211910","J4623020","J4742030","J5102210","J5412110")
  surfaceBV <- c(339,416,164,260,59,25,141,280,117,90,205,165,578,300,620)
  
  i=1
  for (co in code) {
    load(paste0("data/simulations_GR/CemaNeigeGR6J_",co,"_Cal_Compl_Per.RData")) ## change the path !
    data_sav$Qsim<-data_sav$Qsim*surfaceBV[i]/86.4 #mm/day --> m3/s conversion
    assign(paste0("river",i),data_sav)
    i=i+1
  }
  
  ###assemble and select the data
 # library(tibble)
  for (river in 1:length(pops)) {
    river_now<-get(paste0("river",river))
    river_now[,"river"]<-river
    mydates = as.POSIXlt(river_now$DatesR)
    river_now[,"year"]<-format(as.Date(mydates), "%Y")
    river_now[,"month"]<-format(as.Date(mydates), "%m")
    river_now[,"day"]<-format(as.Date(mydates), "%d")
    river_now[,"doy"]<-mydates$yday+1
    
    
    # assign(paste0("river",river),add_column(get(paste0("river",river)),river))
    # 
    # mydates = as.POSIXlt(get(paste0("river",river))$DatesR)
    # assign(paste0("river",river),add_column(get(paste0("river",river)),format(as.Date(mydates), "%Y")))
    # assign(paste0("river",river),add_column(get(paste0("river",river)),format(as.Date(mydates), "%m")))
    # assign(paste0("river",river),add_column(get(paste0("river",river)),format(as.Date(mydates), "%d")))
    # 
    # assign(paste0("river",river),add_column(get(paste0("river",river)),mydates$yday+1))

    if(river==1) {
      #rivers=get(paste0("river",river))
      rivers=river_now
    } else {
      #rivers<-rbind(rivers,get(paste0("river",river)))
      rivers<-rbind(rivers,river_now)
      
    }
  }
  
  colnames(rivers)<-c("date","flow","Qobs","temp","river","year","month","day","doy")
  rivers$year<-as.integer(rivers$year)
  rivers<-rivers[which(rivers$year>=1970&rivers$year<=2017),]

  rivers<-rivers[-which(rivers$month=="02" & rivers$day=="29"),]
  
  ###compute the module
  #mean by month
  monthly_mean <- aggregate(flow~year+month+river, rivers, mean)
  monthly_mean$month<-as.numeric(monthly_mean$month)
  #ponderate mean by year
  yearly_mean_ponderated <- array(,dim=c(length(min(monthly_mean$year):max(monthly_mean$year)),length(pops)+1))
  yearly_mean_ponderated[,1]<-c(min(monthly_mean$year):max(monthly_mean$year))
  for (river in 1:length(pops)) {
    i=1
    for (year in min(monthly_mean$year):max(monthly_mean$year)) {
      yearly_mean_ponderated[i,river+1]<-(monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==01 & monthly_mean$river==river)]*31
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==02 & monthly_mean$river==river)]*28
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==03 & monthly_mean$river==river)]*31
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==04 & monthly_mean$river==river)]*30
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==05 & monthly_mean$river==river)]*31
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==06 & monthly_mean$river==river)]*30
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==07 & monthly_mean$river==river)]*31
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==08 & monthly_mean$river==river)]*31
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==09 & monthly_mean$river==river)]*30
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==10 & monthly_mean$river==river)]*31
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==11 & monthly_mean$river==river)]*30
                                          + monthly_mean$flow[which(monthly_mean$year==year & monthly_mean$month==12 & monthly_mean$river==river)]*31)/365
      i=i+1
    }
  }
  
  #final mean
  yearly_mean_ponderated<-as.data.frame(yearly_mean_ponderated)
  colnames(yearly_mean_ponderated)<-c("year",paste0("river",1:length(pops)))
  
  module=NULL
  for(river in 1:length(pops)) {
    module <- c(module,mean(yearly_mean_ponderated[,river+1]))
  }
  

    #compute Twater
  temperatures=NULL
  flow=NULL
    for (pop in 1:npop){
      tair <- rivers$temp[which(rivers$river==pop)]
      IQ <- 1/((module[pop]/10)+ rivers$flow[which(rivers$river==pop)])
      
      min_Tair_year=(aggregate(temp~year,rivers[which(rivers$river==pop),], min, na.rm=T)$temp)-mean(aggregate(temp~year,rivers[which(rivers$river==pop),], min, na.rm=T)$temp)
      max_Tair_year=(aggregate(temp~year,rivers[which(rivers$river==pop),], max, na.rm=T)$temp)-mean(aggregate(temp~year,rivers[which(rivers$river==pop),], max, na.rm=T)$temp)
      min_flow_year=(aggregate(flow~year,rivers[which(rivers$river==pop),], min, na.rm=T)$flow)-mean(aggregate(flow~year,rivers[which(rivers$river==pop),], min, na.rm=T)$flow)
      #max_flow_year=(aggregate(flow~year,rivers[which(rivers$river==pop),], max, na.rm=T)$flow)-mean(aggregate(flow~year,rivers[which(rivers$river==pop),], max, na.rm=T)$flow)
      
      min_Tair=max_Tair=min_flow=max_flow=NULL
      for (i in 1:length(min_Tair_year)) {
        min_Tair=c(min_Tair, rep(min_Tair_year[i],365))
        max_Tair=c(max_Tair, rep(max_Tair_year[i],365))
        min_flow=c(min_flow, rep(min_flow_year[i],365))
       # max_flow=c(max_flow, rep(max_flow_year[i],365))
      }
      
      #min
      mu.theta <- theta1 + (theta2 * min_Tair) + (theta3 * min_flow)
      theta <- mu.theta + rnorm(1,0,sd.theta)
      #max
      mu.alpha <- alpha1 + (alpha2 * max_Tair) #+ (alpha3 * min_flow)
      alpha <- mu.alpha + rnorm(1,0,sd.alpha)
         
      nyear=nInit + nYears + 1
      nday <- 365
      day <- 1:(365 * nyear)
      
      #compute water temperature
      x=NULL
      x[1] <- sum(a[1] * tair[1]) + b * IQ[1] + c
      x[2] <- sum(c(a[1] * tair[2], a[2] * tair[1])) + b * IQ[2] + c
      x[3] <- sum(c(a[1] * tair[3], a[2] * tair[2], a[3] * tair[1])) + b * IQ[3] + c
      x[4] <- sum(c(a[1] * tair[4], a[2] * tair[3], a[3] * tair[2], a[4] * tair[1])) + b * IQ[4] + c
      x[5] <- sum(c(a[1] * tair[5], a[2] * tair[4], a[3] * tair[3], a[4] * tair[2], a[5] * tair[1])) + b * IQ[5] + c
      x[6] <- sum(c(a[1] * tair[6], a[2] * tair[5], a[3] * tair[4], a[4] * tair[3], a[5] * tair[2], a[6] * tair[1])) + b * IQ[6] + c
      x[7] <- sum(c(a[1] * tair[7], a[2] * tair[6], a[3] * tair[5], a[4] * tair[4], a[5] * tair[3], a[6] * tair[2], a[7] * tair[1])) + b * IQ[7] + c
      x[8] <- sum(c(a[1] * tair[8], a[2] * tair[7], a[3] * tair[6], a[4] * tair[5], a[5] * tair[4], a[6] * tair[3], a[7] * tair[2], a[8] * tair[1])) + b * IQ[8] + c
      x[9] <- sum(c(a[1] * tair[9], a[2] * tair[8], a[3] * tair[7], a[4] * tair[6], a[5] * tair[5], a[6] * tair[4], a[7] * tair[3], a[8] * tair[2], a[9] * tair[1])) + b * IQ[9] + c
      x[10] <- sum(c(a[1] * tair[10], a[2] * tair[9], a[3] * tair[8], a[4] * tair[7], a[5] * tair[6], a[6] * tair[5], a[7] * tair[4], a[8] * tair[3], a[9] * tair[2], a[10] * tair[1])) + b * IQ[10] + c
      x[11] <- sum(c(a[1] * tair[11], a[2] * tair[10], a[3] * tair[9], a[4] * tair[8], a[5] * tair[7], a[6] * tair[6], a[7] * tair[5], a[8] * tair[4], a[9] * tair[3], a[10] * tair[2], a[11] * tair[1])) + b * IQ[11] + c
      x[12] <- sum(c(a[1] * tair[12], a[2] * tair[11], a[3] * tair[10], a[4] * tair[9], a[5] * tair[8], a[6] * tair[7], a[7] * tair[6], a[8] * tair[5], a[9] * tair[4], a[10] * tair[3], a[11] * tair[2], a[12] * tair[1])) + b * IQ[12] + c
      x[13] <- sum(c(a[1] * tair[13], a[2] * tair[12], a[3] * tair[11], a[4] * tair[10], a[5] * tair[9], a[6] * tair[8], a[7] * tair[7], a[8] * tair[6], a[9] * tair[5], a[10] * tair[4], a[11] * tair[3], a[12] * tair[2], a[13] * tair[1])) + b * IQ[13] + c
      x[14] <- sum(c(a[1] * tair[14], a[2] * tair[13], a[3] * tair[12], a[4] * tair[11], a[5] * tair[10], a[6] * tair[9], a[7] * tair[8], a[8] * tair[7], a[9] * tair[6], a[10] * tair[5], a[11] * tair[4], a[12] * tair[3], a[13] * tair[2], a[14] * tair[1])) + b * IQ[14] + c
      x[15] <- sum(c(a[1] * tair[15], a[2] * tair[14], a[3] * tair[13], a[4] * tair[12], a[5] * tair[11], a[6] * tair[10], a[7] * tair[9], a[8] * tair[8], a[9] * tair[7], a[10] * tair[6], a[11] * tair[5], a[12] * tair[4], a[13] * tair[3], a[14] * tair[2], a[15] * tair[1])) + b * IQ[15] + c

      muT=tmpT=NULL
      for (t in 1:15) {
        muT[t] <- theta[t] + ((alpha[t] - theta[t])/(1 + exp(-x[t])))
        tmpT[t] <- muT[t] + rnorm(1,0,sigmaTwater)
      }
      
      
      for (t in 16:length(day)) {
        x[t] <- sum(c(a[1] * tair[t], a[2] * tair[t-1], a[3] * tair[t-2], a[4] * tair[t-3], a[5] * tair[t-4], a[6] * tair[t-5], a[7] * tair[t-6], a[8] * tair[t-7], a[9] * tair[t-8], a[10] * tair[t-9], a[11] * tair[t-10], a[12] * tair[t-11], a[13] * tair[t-12], a[14] * tair[t-13], a[15] * tair[t-14], a[16] * tair[t-15])) + b * IQ[t] + c
        muT[t] <- theta[t] + ((alpha[t] - theta[t])/(1 + exp(-x[t])))
        tmpT[t] <- muT[t] + rnorm(1,0,sigmaTwater)
      }
     
      #tmpT<-temp.simulated
      tmpT[tmpT < 0] <- 0.01
      #tmpT <- tmpT+ 0.5
      # tmpT[tmpT < 5] <- 5
      # tmpT[tmpT > 19.5] <- 19.5
      
      #plot(tmpT[(365*30):(365*40)], type="l")
      
      #compute logrelflow
      #  tmplogrelflow <- log(rivers$flow[which(rivers$river==pop)]/module[pop])
      #  tmplogrelflow <- tmplogrelflow -  mean(tmplogrelflow) #correction due to values not centered on 0
      # # 
      # #tmplogrelflow<-debit.simulated
      #  #tmplogrelflow[tmplogrelflow > 3] <- 3
      #  tmplogrelflow[tmplogrelflow > 1.8] <- 1.8
      #  
      # #tmplogrelflow[tmplogrelflow < -3] <- -3
      # tmplogrelflow[tmplogrelflow < -2.3] <- -2.3
      
      tmpflow <- rivers$flow[which(rivers$river==pop)]
      
      
      temperatures <- cbind(temperatures,tmpT)    
      flow <- cbind(flow,tmpflow)
      
    } # end loop pop
    
    #return(list(mTair=mTair,ampF=ampF,temperatures = temperatures, logrelflow = logrelflow))
    
  #} #end function river_climate_model_multi
  
 # env <- river_climate_model_multi(npop, nInit, nYears, tempCC, ampCC) 
  env <- list(temperatures=temperatures, flow=flow, module = module)
  
} #end if Obs=TRUE


# 
# ### visualisation
# plot(temperatures[1:(365*2),1], type='l')
# plot(logrelflow[1:(365*2),1], type='l')
# #nbr of time sup or inf compared to critic flow (once data are centered around 0)
# cfi=log(0.1)
# cfs=log(7)
# 
# nb_critinf=NULL
# nb_critsup=NULL
# for (river in 1:length(pops)) {
#   nb_critinf <- c(nb_critinf , length(logrelflow[,river][which(logrelflow[,river]<cfi)]))
#   nb_critsup <- c(nb_critsup , length(logrelflow[,river][which(logrelflow[,river]>cfs)]))
# }
# 
# crit <- cbind(pops,nb_critinf, nb_critsup)
# colnames(crit) <- c("Pop","Nb crit inf", "Nb crit sup")
# 
# mean_temp=NULL
# #mean temperatures
# for (river in 1:length(pops)) {
#   mean_temp <- c(mean_temp , mean(temperatures[,river]))
# }
# temp<-cbind(pops,mean_temp)
# #save(river_climate, file="data/river_climate.Rdata")
# 
# ## PLOT
# png(filename = "data/Climate.png")
# par(mfrow=c(2,2))
# 
# ny = 50
# x0 <- 1 #365*40
# x<- (x0-1) + 365*ny
# 
# ## TEMPERATURES
# plot(NULL, xlim=c(0,x),ylim=range(env$mT[x0:x]),xlab="Years",ylab="Mean Water temperature",xaxt='n')
# for (pop in 1:npop){
#   lines(env$mT[x0:x] ,type='l', col=pop)
#   #lines(mT[x0:x] ,type='l', col=pop)
# }
# mtext(seq(0,ny, 5), side = 1, line = 1, outer = FALSE, at = seq(0,ny, 5)*365)
# 
# 
# plot(NULL, xlim=c(0,x),ylim=range(env$temperatures[x0:x,]),xlab="Days",ylab="Water temperature")
# for (pop in 1:npop){
#   lines(env$temperatures[x0:x,pop], ,type='l', col=pop)
# }
# 
# ## FLOW
# plot(NULL, xlim=c(0,x),ylim=range(env$ampF[x0:x]),xlab="Years",ylab=" Mean amplitude Flow (m3/s)",xaxt='n')
# for (pop in 1:npop){
#   lines(env$ampF[x0:x] ,type='l', col=pop)
#   #lines(mT[x0:x] ,type='l', col=pop)
# }
# mtext(seq(0,ny, 5), side = 1, line = 1, outer = FALSE, at = seq(0,ny, 5)*365)
# 
# plot(NULL, xlim=c(0,x),ylim=range(env$logrelflow[x0:x,]),xlab="Days",ylab="Flow (m3/s)")
# for (pop in 1:npop){
#   lines(env$logrelflow[x0:x,pop], ,type='l', col=pop)
# }
# dev.off()
