#set.seed(666) # R‘s random number generator. Use the set.seed function when running simulations to ensure all results, figures, etc are reproducible.

#_________________ GENERAL PARAMETERS _________________#
#pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
#pop = 15 # population to simulate
#npop = length(pops)
nYears = 40 #for data observed: 37
nInit = 10
rPROP = .25


#_________________  DATA POPULATIONS _________________#
## 1.1 Population characteristics
dat <- read.csv2("data/dataPop.csv", header=TRUE, stringsAsFactors = FALSE, comment.char = "#")
#attach(dat)
#pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
#dat <- dat[which(dat$Population %in% pops),] # remove Couesnon

dat <- dat[-which(dat$Population == "Couesnon"),] # remove Couesnon
#dat <- dat[which(dat$Population == "Scorff"|dat$Population == "Yar"),] # remove Couesnon

#edit al - 20/03/21 - for all pop the same = Scorff
# test<-dat[which(dat$Population == "Scorff"),]
# for (i in 1:14) {
#   test <- rbind(test,dat[which(dat$Population == "Scorff"),])
# }
# dat<-test
#

pops <- dat$Population; npop <- length(pops)
M=as.numeric(dat$M)

#scn SRR variable
# load("/media/hdd/alamarins/metapop_amaia/7.third-diversity-design/data/pops_contrasted_Diversity_TRUE.RData")
# pops_contrasted[which(pops_contrasted=="Ellee")] <- "Laita"
# 
# M <- rep(1,15)
# for (pop in pops_contrasted[1:5]) { #+2° for 5 first pops
#   id = which(pops==pop)
#   M[id] <- 2.6#1.2
# }
# for (pop in pops_contrasted[6:10]) { #-2° for 5 last pops
#   id = which(pops==pop)
#   M[id] <- 0.4#0.8
# }

scenarioSize = as.numeric(args[8])
if (scenarioSize == 1) { #equal size btw pop
  Area=rep(sum(dat$Area)/npop, npop)
}
if (scenarioSize == 2) { #"real" size
  Area=as.numeric(dat$Area)
}

Rmax=as.numeric(dat$Rmax)
Alpha=as.numeric(dat$alpha)
Type=dat$Type



#_________________ SCENARIOS _________________#



############# 1. CONNECTIVITY #############

# 1.1 Parameters for Laplace kernel
mu=0
beta=29.5 # so that >80% of migrants disperse into the first 50km

h = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5) # Philopatry (homing) rates #edit al - 22/03/21 5% and 15%
#scenarioConnect=7 #scenario 1 for h=1.00, scenario 2 for h=0.95, scenario 5 pour h=0.80
scenarioConnect = as.numeric(args[1]) #scenario 1 for h=1.00, scenario 2 for h=0.95, scenario 3 pour h=0.80

# 1.2 Connectivity matrix
source("code/Connectivity_matrix.R")
pstray = connect_kernel
# if (npop < 2){
#   pstray = matrix(1,npop,npop)
# } else {
#   # première matrice h=1, deuxième matrice h=0.942, troisième matrice h=0.80
#   #load("data/Matrices_Laplace_AireLog.RData")
#   tmp = connect[[scenarioConnect]]
#   tmp = tmp[,which(colnames(tmp) %in% toupper(pops))]
#   pstray = tmp[which(rownames(tmp) %in% toupper(pops)),]
# }



############ 2. ENVIRONMENT ############

## 2.1 Environmental conditions (discharge, water temperature)


#1. CHOOSE WHICH TYPES OF ENVIRONMENTAL DATA YOU WANT

Obs=FALSE #here you can choose to use observed environmental data (flow in mm/day and air temperature) or to simulate them

#2. if you want to simulate them:
Diversity_env=FALSE #TRUE #choose if you want different regimes between rivers (TRUE) or same regime as Scorff (TRUE)

## Climate change scenarios
scenarioEnvi=1 #scenario 1 pour absence de CC, scenario 2 pour CC
env=list(c(0,1,1), c(3,1.25,0.75))
tempCC=env[[scenarioEnvi]][1]
ampCC=env[[scenarioEnvi]][2]
seaCC=env[[scenarioEnvi]][3]

rhoT<-0 # correlations temperatures fluctuations
corMat <- array(rhoT,dim=c(npop,npop)); diag(corMat)<-1
rhoF<-0 # correlations flow fluctuations
#corMat <- array(rhoF,dim=c(npop,npop)); diag(corMat)<-1


#3. if you want to use the same conditions for several scenarios
# 3.1. generate them from the generate_env_conditions.R code for N simulations
# 3.2. if already generated, load them from the data folder: env_new_false.Rdata & env_new_true.Rdata
if (file.exists(paste0("data/environmental_conditions_Diversity_",Diversity_env,"_",last,"simul.RData"))) {
  load(paste0("data/environmental_conditions_Diversity_",Diversity_env,"_",last,"simul.RData"))
} else {
  source("code/generate_env_conditions.R")
}

#load("data/environmental_conditions_Diversity_TRUE_100simul.RData")

#4.if want to generate new environmental conditions for each scenario, use the river_climate_multi_Tair script in the simulation loop in metaIbasam.R





################ 3. FISHERIES ###############@
fish.state=TRUE # TRUE If fishing applied
fish.stage=TRUE # fishing on life stages (1SW/MSW) if TRUE, on Sizes ("small","med","big") otherwise

## 3.1. Fishing rate

#by stage
frates_vec1SW =c(0, 0.05, 0.07, 0.1)
frates_vecMSW =c(0, 0.1, 0.15, 0.2)

#equal between stages
#frates_vec =c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
scenarioFishingRate = as.numeric(args[3])


#frates=c(.0,.15,.30) # Fishing rates
if ( fish.stage==TRUE) {
  frates = c(0, frates_vec1SW[scenarioFishingRate], 0, frates_vecMSW[scenarioFishingRate]) # (1SW_init, 1SW, MSW_init, MSW)
  
  fratesSource = frates # fishing rates by stage (1SW vs MSW) for source populations (13/04/2021 : ajout du taux d'exploitation initial et au bout de 10ans) 
  fratesSink = frates # fishing rates by stage (1SW vs MSW) for sink populations
  fratesNeutral = frates # fishing rates by stage (1SW vs MSW) for neutral populations
}

if ( fish.stage==FALSE) {
  frates = c(0.07, 0.10, 0.15) # fishing rates by size
}


## 3.2. Fishing scenarios
# 0: control /
# 1: no fishing on sink populations
# 2: no fishing on source populations
# 3: no fishing on sink populations & stronger fishing pressure on sources
# 4: no fishing on source populations & stronger fishing pressure on sinks
scenarioFishing = as.numeric(args[2])

# fishing rates (frates) are provided in parIbasam.R file
tmp <- matrix(0, nrow = npop, ncol = length(fratesSink))
if (scenarioFishing==0){ # fishing all
  for (i in 1:npop) {
    if (dat$Type[i]=="sink") tmp[i,] <- fratesSink
    if (dat$Type[i]=="neutral") tmp[i,] <- fratesNeutral
    if (dat$Type[i]=="source") tmp[i,] <- fratesSource
  }}

if (scenarioFishing==1){ # fishing sources/neutral
  for (i in 1:npop) {
    if (dat$Type[i]=="sink") tmp[i,] <- 0
    if (dat$Type[i]=="neutral") tmp[i,] <- fratesNeutral
    if (dat$Type[i]=="source") tmp[i,] <- fratesSource
  }}

if (scenarioFishing==2){ # fishing sinks/neutral
  for (i in 1:npop) {
    if (dat$Type[i]=="sink") tmp[i,] <- fratesSink
    if (dat$Type[i]=="neutral") tmp[i,] <- fratesNeutral
    if (dat$Type[i]=="source") tmp[i,] <- 0
  }}

# if (scenarioFishing==3){ # fishing sources + neutral
#   for (i in 1:npop) {
#     if (dat$Type[i]=="sink") tmp[i,] <- 0
#     if (dat$Type[i]=="neutral") tmp[i,] <- frates * 1.5 # increase fishing rates by 2 (compensation)
#     if (dat$Type[i]=="source") tmp[i,] <- frates * 1.5
#   }}
# 
# if (scenarioFishing==4){ # fishing sinks + neutral
#   for (i in 1:npop) {
#     if (dat$Type[i]=="sink") tmp[i,] <- frates * 1.5
#     if (dat$Type[i]=="neutral") tmp[i,] <- frates * 1.5
#     if (dat$Type[i]=="source") tmp[i,] <- 0
#   }}
fishing.rates <- tmp

# 3 scenarios (iSIMUL) for STAT:
#grilse: 0.15 / MSW: 0.15 => "Control"
#grilse: 0.3 / MSW: 0 => "1SW"
#grilse: 0. / MSW: 0.3 => "MSW"

# 5 scenarios (iSIMUL) for TAILLE:
# small: 0.15 / med: 0.15 / big: 0.15 => "Control"
# small: 0.6 / med: 0 / big: 0 => "Small"
# small: 0 / med: 0.3 / big: 0 = "Med"
# small: 0 / med: 0 / big: 0.6 => "Big"
# small: 0.3 / med: 0 / big: 0.3 => "Big small"




################ 4. GENETIC DIVERSITY ###############@

## GENETIC CHANGES AT INITIALIZATION - al 12/02/2021
#propG=rep(0,npop)
#propTO=rep(0,npop)

diversity_prop = c(0,0.1,0.15)

scenarioDiversityG = as.numeric(args[4]) #scn of genetic change at initialization
scenarioDiversityTO = as.numeric(args[5]) #scn of growth survival trade-off parameter sigRIV change


if (scenarioDiversityG == 1) { #no diversity - all pops at +0.15
  propG <- seq(from = diversity_prop[scenarioDiversityG], to = -diversity_prop[scenarioDiversityG], length.out = 15) + 0.15
} else {
  propG <- seq(from = diversity_prop[scenarioDiversityG], to = -diversity_prop[scenarioDiversityG], length.out = 15) #gradient of diversity
}

propTO <- seq(from = diversity_prop[scenarioDiversityTO], to = -diversity_prop[scenarioDiversityTO], length.out = 15)

scenarioDiversityStruct = as.numeric(args[6])

if (scenarioDiversityStruct == 2) { #random diversity -> shuffle
  #shuffle propG
  #random_propG <- sample(propG)
  #save(random_propG, file=paste0("data/random_propG_", diversity_prop[scenarioDiversityG],".RData"))
  #random_propTO <- sample(propTO)
  #save(random_propTO, file=paste0("data/random_propTO_", diversity_prop[scenarioDiversityTO],".RData"))
  
  #load randomized propG
  load(paste0("data/random_propG_",diversity_prop[scenarioDiversityG],".RData"))
  propG <- random_propG
  #load randomized propTO
  load(paste0("data/random_propTO_",diversity_prop[scenarioDiversityTO],".RData"))
  propTO <- random_propTO
}

# if (scenarioDiversityStruct == 1) { #structured diversity
#   for (pop in 1:5) {
#     propG[pop] <- diversity_prop[scenarioDiversityG]
#     propTO[pop] <- diversity_prop[scenarioDiversityTO]
#   }
#   for (pop in 11:15) {
#     propG[pop] <- -diversity_prop[scenarioDiversityG]
#     propTO[pop] <- -diversity_prop[scenarioDiversityTO]
#   }
# }
  
# if (scenarioDiversityStruct == 2) { #random diversity
#   #load pops contrasted
#   load("data/pops_contrasted_Diversity_TRUE.RData")
#   pops_contrasted[which(pops_contrasted=="Ellee")] <- "Laita"
#   pops[which(pops=="Ellee")] <- "Laita"
#   
#   for (pop in pops_contrasted[1:5]) { #+X% for 5 first pops
#     id = which(pops==pop)
#     propG[id] <- diversity_prop[scenarioDiversityG]
#     propTO[id] <- diversity_prop[scenarioDiversityTO]
#   }
#   for (pop in pops_contrasted[6:10]) { #-X% for 5 last pops
#     id = which(pops==pop)
#     propG[id] <- -diversity_prop[scenarioDiversityG]
#     propTO[id] <- -diversity_prop[scenarioDiversityTO]
#   }
# }



