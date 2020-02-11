install.packages()

### Community simulation

rm(list=ls())
library(igraph)

library(magrittr)
library(purrr)
library(readr)
library(matrixcalc)
library(plot.matrix)
library(devtools)
library(vegan)
#install_github('nathanvan/parallelsugar')
library(parallelsugar)


workdir= "C:/Users/admbotella/Documents/pCloud local/boulot/data/Simu_Networks/"


simulate_community <- function(
  env = runif(50, 0, 100),
  niche_optima  = seq(2, 98, length.out=50), 
  niche_breadth = 20,
  comp_inter = NA, 
  fac_inter = NA, 
  beta_env = 1,
  beta_comp = 5, 
  beta_fac = 0, 
  beta_abun = 0, 
  years = 100, 
  K = 40,
  competition = "facilitation", intra_sp_com  = 0, epochs=T
  
) {
  #code for one site
  sim_com <- function(
    env, niche_breadth, niche_optima, comp_inter, fac_inter, beta_env,
    beta_comp, beta_fac, beta_abun, years, K, competition,intra_sp_com, epochs
  ) {
    n_sp = length(niche_optima)
    
    # Define adjacency matrices for negative ("Competition") and positive ("Facilitation") interactions
    # Interactions doesn't have to be symmetric, i.e. assymetric adjacency matrix allowed.
    if (length(comp_inter) == 1) comp_inter = matrix(comp_inter, n_sp, n_sp)
    if (length(fac_inter)  == 1) fac_inter  = matrix(fac_inter, n_sp, n_sp)
    species_niche_overlap_sym <- comp_inter #Negative interaction matrix
    species_fac_sym <- fac_inter #Positive interaction matrix
    
    
    log_p_env <- sapply(niche_optima, dnorm, mean = env, sd = niche_breadth, log = TRUE)
    log_p_env <- log_p_env  - log(dnorm(0) / 10)
    
    community <- factor(
      x      = sample(seq_along(niche_optima), K, replace = TRUE),
      levels = seq_len(n_sp)
    )
    
    abund <- table(community)
    hist<- matrix(NA, nrow=years,ncol=n_sp)
    
    
    for (j in seq_len(years)) {
      for (k in seq_len(K)) {
        PositiveInteraction <- 1 - colSums(species_fac_sym[community, ]) / K
        NegativeInteraction <- 1 - colSums(species_niche_overlap_sym[community, ]) / K
        
        # Probabiliy of each species to enter the community
        if (competition == "facilitation") {
          p_all <- exp(
            beta_env * log_p_env - beta_fac * log(PositiveInteraction) + 
              log(1 + beta_abun * abund)
          )
        } else if(competition=="all"){
          p_all <- exp(
            beta_env * log_p_env + beta_comp * log(NegativeInteraction)- beta_fac * log(PositiveInteraction) +
              log(1 + beta_abun * abund)
          )
        }else if(competition=="competition"){
          p_all <- exp(
            beta_env * log_p_env + beta_comp * log(NegativeInteraction) +
              log(1 + beta_abun * abund)
          )
        }
        
        
        p_all <- ifelse(is.na(p_all), min(p_all, na.rm = TRUE), p_all)
        if (all(is.na(p_all)) || identical(min(p_all), max(p_all))) p_all = NULL
        if (any(is.infinite(p_all))) {
          community[sample(K, 1)] <- sample(seq_len(n_sp)[p_all == Inf], 1)
        } else {
          community[sample(K, 1)] <- sample(n_sp, 1, prob = p_all) #K times an individual is randomly selected
        }                                                          #to be replaced by an individual selected by p_all
        #p_all updates over K
        abund <- table(community)
        
      }
      
      hist[j,]<-abund
    }
    if(!epochs) {
      as.integer(abund) > 0 #if P/A
      #abund #if counts
    }else{
      #if all the histogram
      as.data.frame(hist)
    }
  }
  
  ans <- parallelsugar::mclapply(
    env, sim_com, niche_breadth, niche_optima, comp_inter, fac_inter,
    beta_env, beta_comp, beta_fac, beta_abun, years, K, competition,
    intra_sp_com,epochs
    #,mc.cores = detectCores()
  )
  
  if(!epochs){
    #if just last iteration
    ans <- do.call(rbind, ans)
    ans <- cbind(ans, env)
    sp_labs <- paste0(
      "sp_", gsub(" ", 0, format(seq_along(niche_optima), width = 2))
    )
    colnames(ans) <- c(sp_labs, "env")
    as.data.frame(ans)
  }else{
    #if all epochs
    sp_labs <- paste0(
      "sp_", gsub(" ", 0, format(seq_along(niche_optima), width = 2))
    )
    for(i in 1:length(ans)) colnames(ans[[i]])<-sp_labs
    #lapply(ans, function(y) colnames(y) <- sp_labs)
    list(hist=ans, env=env)
  }
}

#just a test to see if the function works
#test.coms <- simulate_community()

### Function to create the trophic interaction matrix

# It follows the niche model of Williams and Martinez (2000)
# User must specify
# @S : number of species
# @minR : Minimum trophic range. The connectance is all the bigger as it is highand it favors loops
makeTrophicInteractionMatrix=function(S,minR = 0.2){
  n_is = runif(S,0,1)
  r_is = runif(S,minR,1)
  r_is = r_is * n_is
  c_is = runif(S,0,1)
  c_is = c_is * n_is
  
  A = matrix(0,S,S)
  for(i in 1:S){
    # (i,j) equals -1 if on site presence of i bother the presence of j 
    # Thus the -1 on the line i indicate all the columns that are eaten by it 
    A[i,] = A[i,] -1*as.numeric(n_is>=(c_is[i]-r_is[i]/2) & n_is<=(c_is[i]+r_is[i]/2))
    # (i,j) equals +1 if on site presence of i favors the presence of j
    # Thus the +1 on the line i indicate all the columns that eat i
    A[i,] = A[i,] + as.numeric((c_is - r_is/2)<=n_is[i] & (c_is + r_is/2)>=n_is[i])
  }
  return(A)
}


### Creation of the matrices
#We create the adjacency matrices with different parameters. In particular, we vary the type of interactions (facilitation or competition), the percentage of interactions that we sample among all the ones that follow the criteria and the number of species. We plot the so-built adjacency matrix.

# Sp is 50 or 100
# We take niche_breath=20. This means that with 50 species, there are 345 possible interactions (each species can compete/facilitate with 10 species to the left and 10 to the right). With 100 species, there are 1449 species, and each species can compete/facilitate with 20 species to the left and 20 to the right.
# We take Nint = 10, 20, 50, 100

### Simulation 

Ss = c(5,10,20,30)
minRs = c(.1,.2,.3,.4)

df = expand.grid(S=Ss,minR=minRs)


sim_names= paste('S',df$S,'minR',df$minR,sep="")

nruns=length(sim_names)
set.seed(1123)


nruns=length(sim_names)
sim_params =list(S= df$S,minR= df$minR) 

adj_mats = sim_params %>% pmap(makeTrophicInteractionMatrix)
adj_mats=  adj_mats %>% set_names(sim_names) 


### Plot the Interaction networks
setwd(workdir)
if(!"AdjMatrices"%in%list.files()){dir.create("AdjMatrices")}
setwd(paste(workdir,'AdjMatrices/',sep=""))
if(!"Plots"%in%list.files()){dir.create("Plots")}
lapply(sim_names, function(x){
  png(file=paste0(paste(workdir,"AdjMatrices/Plots/",x,sep=""),".png"),height=500,width=500)
  plot.igraph(graph_from_adjacency_matrix(adj_mats[[x]]))
  dev.off()
})

saveRDS(adj_mats,file=paste(workdir,"AdjMatrices/adjmatrices.rds",sep=""))


### Creating communities
#We create the communities following the parameters above. For each simulation, we take the interaction network created above.

setwd(workdir)
if(!"Communities"%in%list.files()){dir.create("Communities")}
outputdir=paste(workdir,"Communities/",sep="")



As= readRDS(file=paste(workdir,"AdjMatrices/adjmatrices.rds",sep=""))

lcomp= lapply(As , function(A){
  Aneg = A
  Aneg[A==-1] = 0
  return(Aneg)
})

lfac = lapply(As , function(A){
  Apos = -A
  Apos[A==-1] = 0
  return(Apos)
})


sim_params =list(
  niche_optima = rep( map( Ss , function(x) seq(2,98,length.out= x)) , round(nruns/length(Ss)) ),
  comp_inter   = lcomp,
  fac_inter    = lfac,
  beta_comp    = rep(list(10),each=nruns),
  beta_fac     = rep(list(5),each=nruns),
  beta_env     = rep(list(1), each = nruns),
  beta_abun    = rep(list(5), each = nruns),
  K            = rep(list(40), each = nruns),
  competition  = rep(list("all"),each=nruns),
  intra_sp_com = rep(list(0), nruns),
  epochs       = rep(T,nruns)
) 
sim_data = sim_params %>% pmap(simulate_community) 
sim_data=sim_data %>% set_names(sim_names) 
saveRDS(sim_data,file=paste(outputdir,"sim_data.rds",sep=""))

### Convergence validation
#Each boxplot represent the distribution of the shannon index over all communities at one time step. All the plots show convergence for any kind of parameter when we take 100 time steps.


alpha_div=function(tab){ ##tab with rows=epochs and columns=species
  a=vegan::diversity(tab,"shannon")
}

### nruns = nombre de simulations
### The following plots evolution of alpha-diversity (given by Shannon index applied to species counts) over epochs
### Each iteration corresponds to one simulation model

for (r in 1:nruns){
  hist=sim_data[[r]]$hist  ##contains the data for each epoch rows=epochs and columns=species
  df=do.call(rbind,lapply(hist,alpha_div)) #sites x years matrix
  png(file=paste0(paste(outputdir,"ConvPlots/",sim_names[[r]],sep=""),".png"),height=800,width=1400)
  boxplot.matrix(df,main=paste("Alpha diversity evolution over epochs: ",sim_names[r]),xlab="Simulation      epochs",ylab="Average alpha diversity") 
  dev.off()
}


### Save results
lapply(sim_names, function(x){
  hist=sim_data[[x]]$hist
  env=data.frame(sim_data[[x]]$env)
  #env_poly<-poly(env,2) #if we want to give the poly of env
  occur=do.call(rbind,lapply(hist,function(y) y[nrow(y),]))  ###Keeping only the last community composition
  occur=as.data.frame(do.call(cbind,lapply(occur,function(y) as.integer(y>0)))) # Set to PA
  occur[,"env"]<-env
  occur = occur[order(occur$env),,drop=F]
  write.csv2(occur,file=paste(outputdir,paste0(x,".csv"),sep="/"))
})
