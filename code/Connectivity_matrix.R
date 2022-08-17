##### CONNECTIVITY MATRIX #####

#_________________ PARAMETERS _________________#
#pops <- c("Leff", "Trieux", "Jaudy", "Leguer", "Yar", "Douron", "Penze", "Elorn", "Aulne", "Goyen", "Odet", "Aven", "Laita", "Scorff", "Blavet")
#pop = 15 # population to simulate
#npop = length(pops)
#source("parIbasam.R")


##### I. DISTANCE ####
distance=matrix(data=NA, ncol=npop, nrow=npop)
colnames(distance)=rownames(distance)=pops

scenarioDist = as.numeric(args[7])

if (scenarioDist == 1) {
  tot_length <- dat$Distance[15] - dat$Distance[1]
  dist<-NULL
  for (i in c(0:14)) {
    dist <- c(dist, tot_length/(npop-1) * i)
  }
  for (i in 1:npop){
    for (j in 1:npop){
      distance[i,j]=abs(dist[j]-dist[i])
    }}
}
if (scenarioDist == 2) {
  dist <- dat$Distance  #edit al - 20/03/21
  #dist <- seq(from=600, length.out = 15, by=40) #edit al - 20/03/21 - all pop same distance (40km)
  
  #diff.dist <- c(0,diff(dat$Distance))
  for (i in 1:npop){
    for (j in 1:npop){
      distance[i,j]=abs(dist[j]-dist[i])
    }}
}


# dist_normees=matrix(data=NA, ncol=16, nrow=16)
# colnames(dist_normees)=colnames(distance)
# rownames(dist_normees)=rownames(distance)
# for (i in 1:16){
#   for (j in 1:16){
#     dist_normees[i,j]=(distance[i,j]-dist_moy_sigma[i,"Moyenne"])/dist_moy_sigma[i,"Ecart"]
#   }
# }
#write.table(dist_normees,file="Distances_normees.txt")



##### II. AREA #####
area_log=log10(Area)
#comprendre ce log10
#log10(20) doit �tre deux fois inf�rieur � log10(40)
#log10(40) doit �tre deux fois inf�rieur � log10(20)
#log10(c(200, 400, 800))
#m�me diff�rence entre log10(200) et log10(400) et entre log10(400) et log10(800)

ratio_area=matrix(, nrow=1, ncol= npop); colnames(ratio_area)=pops
for (i in 1:npop){
  ratio_area[1,i]=area_log[i]/(sum(area_log))} 


##### III. DISPERSAL KERNEL #####
lap = function(mu, beta, distance) {(1/(2*beta))*exp(-(distance-mu)/beta)} # Laplace
#dist_Lap=seq(from=4, to= 370, by=23) 
#ressemble au vecteur distance r�el que j'ai
#16 distances avec 50km pour avoir p(d<=50km)>80%(1-h)




##### IV. CONNECTIVITY MATRIX #####

connect_kernel=matrix(, nrow=npop, ncol=npop)
rownames(connect_kernel)=colnames(connect_kernel)=pops
#diag(connect_kernel)=h[scenarioConnect]
connect=list()

for (i in 1:npop){
  attrac_d=lap(mu, beta, distance[i,]) #attractivity score of sites based on distance
  attrac_d[i]<-0 # ignore population of origin
  attrac_ad=attrac_d*ratio_area # influence of area on attractivity
  proba_attract_ad=attrac_ad/sum(attrac_ad) #proba of attractivity du site fonction distance at aire
  sum(proba_attract_ad)
  
  proba_dispersion=proba_attract_ad*(1-h[scenarioConnect]) #proba de dispersion de Couesnon vers autre pop
  proba_dispersion[i] <- h[scenarioConnect]
  connect_kernel[i,]=proba_dispersion #je remplis ma matrice pour Couesnon en tant que pop origine
  test <- rowSums(connect_kernel) #sum of probabilities / need to be 1
  if (!(any(test==1))) {
    print ("something wrong with dispersal kernel")
    connect_kernel <- NULL
  }
  #colSums(connect_kernel)
}
#connect[[paste0(scenarioConnect)]] <- connect_kernel
#save(connect_kernel, file="data/connectivity_matrix.RData")
