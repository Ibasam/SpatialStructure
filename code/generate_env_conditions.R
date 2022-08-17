############### TO GENERATE ENV CONIDITIONS FOR SEVERAL SIMULATIONS AND SCENARIOS


#_________________ PARAMETERS _________________#
#source("parIbasam.R") #definition of number of years, CC scenario, covariation between pops

nSIMUL=as.numeric(args[10])#100

#_________________________________________________#


if (Diversity_env == FALSE) { # FOR SCN NO DIVERSITY
  env_new_false<-list()
  for (simul in 1:nSIMUL) {
    source("code/river_climate_multi_Tair.R")
    
    env_new_false[[simul]]<-list(temperatures=env$temperatures, flow=env$flow, module=env$module)
    
  }
  save(env_new_false,file=paste0("data/environmental_conditions_Diversity_FALSE_",nSIMUL,"simul.RData"))
  
} else { #FOR SCN REAL DIVERSITY
  env_new_true<-list()
  for (simul in 1:nSIMUL) {
    source("code/river_climate_multi_Tair.R")
    #env$temperatures[,14]<-env_new_false[[simul]]$temperatures[,14]
    #env$flow[,14]<-env_new_false[[simul]]$flow[,14]
    env_new_true[[simul]]<-list(temperatures=env$temperatures, flow=env$flow, module=module)
  }
  save(env_new_true,file=paste0("data/environmental_conditions_Diversity_TRUE_",nSIMUL,"simul.RData"))
  
}


# #_________________________________________________#
# #FOR SCN DIVERSITY CONTRASTED
# #random choice of which pop contrasted (10 over 15, except Scorff)
# #first 5 will be +2°, second 5 will be -2°
# probabilities=c(rep(1,13),0,1)
# pops_contrasted<-sample(pops,10,prob=probabilities)
# save(pops_contrasted,file=paste0("data/pops_contrasted_Diversity_TRUE.RData"))
# Diversity_env=TRUE #FALSE
# #not generate new chroniques of env conditions to keep the same as the control scn and just change the mean temperature/flow
# 
# ## /!\/!\/!\/!\/!\ change Ellee by Laita in pops_contrasterd /!\/!\ !!!!!!! 
# #(if need to generate new env conditions)
# pops_contrasted[which(pops_contrasted=="Ellee")] <- "Laita"
# 
# ### TEMPERATURE 
# env_new_true<-list()
# for (simul in 1:nSIMUL) {
#   env_new_true[[simul]] <- env_new_false[[simul]] #same flow for all pop and same temp for intermediaire pop
#   for (pop in pops_contrasted[1:5]) { #+2° for 5 first pops
#     id = which(pops==pop)
#     env_new_true[[simul]]$temperatures[,id] <- env_new_true[[simul]]$temperatures[,id] + 2
#   }
#   for (pop in pops_contrasted[6:10]) { #-2° for 5 last pops
#     id = which(pops==pop)
#     env_new_true[[simul]]$temperatures[,id] <- env_new_true[[simul]]$temperatures[,id] - 2
#   }
# }
# save(env_new_true,file=paste0("data/environmental_conditions_Diversity_TRUE_100simul.RData"))
# 
# ### FLOW 
# env_new_true<-list()
# for (simul in 1:nSIMUL) {
#   env_new_true[[simul]] <- env_new_false[[simul]] #same flow for all pop and same temp for intermediaire pop
#   for (pop in pops_contrasted[1:5]) { #+20% for 5 first pops
#     id = which(pops==pop)
#     env_new_true[[simul]]$flow[,id] <- env_new_true[[simul]]$flow[,id] * 1.2
#   }
#   for (pop in pops_contrasted[6:10]) { #-20% for 5 last pops
#     id = which(pops==pop)
#     env_new_true[[simul]]$flow[,id] <- env_new_true[[simul]]$flow[,id] * 0.8
#   }
# }
# save(env_new_true,file=paste0("data/environmental_conditions_Diversity_TRUE_FLOW_30simul.RData"))
# 
