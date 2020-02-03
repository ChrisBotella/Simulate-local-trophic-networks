

### Community simulation

rm(list=ls())
library(magrittr)
library(purrr)
library(readr)
library(matrixcalc)
library(plot.matrix)
library(devtools)
library(vegan)
#install_github('nathanvan/parallelsugar')
library(parallelsugar)
setwd("C:/Users/user/pCloud local/boulot/data/Simu_Networks/")

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
    species_niche_overlap_sym <- comp_inter #competion matrix
    species_fac_sym <- fac_inter #facilitation matrix
    
    
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
        f_comp <- 1 - colSums(species_fac_sym[community, ]) / K
        p_comp <- 1 - colSums(species_niche_overlap_sym[community, ]) / K
        
        # Probabiliy of each species to enter the community
        if (competition == "facilitation") {
          p_all <- exp(
            beta_env * log_p_env - beta_fac * log(f_comp) + 
              log(1 + beta_abun * abund)
          )
        } else {
          p_all <- exp(
            beta_env * log_p_env + beta_comp * log(p_comp) +
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

### Function to create the interactions
# Code to create the interaction network.
# The parameters of such functions are:
# @Sp, the number of species
# @niche_breadth, the standard deviations of the niches
# @type, whether we want to create a facilitation network or a competition network.
# @Nint is the % of all the possible interactions that we sample from the set of all possible interactions.
# For facilitation, two species can interact if their niche are not too close nor too far. That is, if their niche optima mu1,mu2 are s.t. sigma<|mu1-mu2|<2*sigma. We instead suppose that two species can compete if their niches are close, i.e. if |mu1-mu2|<sigma.


CreateInteractionMatrix<-function(Sp=20,Nint=10,niche_breadth = 20,type="facilitation"){

niche_optima  = seq(2, 98, length.out=Sp)

if(!(type=="facilitation" | type=="competition")) stop("Type must be facilitation or competition")
if(Nint<0 | Nint>100) stop("Nint must be a percentage")

if(type=="facilitation") {factors=c(1,2)}else{factors=c(0,1)}

# Two species compete if their niche optima mu1,mu2 are s.t. |mu1-mu2|<sigma
# Two species facilitate if their niche optima mu1,mu2 are s.t. sigma<|mu1-mu2|<2*sigma
all_int_mat <- outer(
        niche_optima,
        niche_optima,
        function(x, y) ifelse(abs(x-y)>factors[1]*niche_breadth & abs(x-y)<factors[2]*niche_breadth,1,0)
      )
all_int_mat[lower.tri(all_int_mat,diag=TRUE)] <- 0
n_poss<-length(which(all_int_mat==1))
interactions<-sample(which(all_int_mat==1),Nint*n_poss/100,replace=F)

sampled_int_asym<-matrix(0,Sp,Sp)
sampled_int_asym[interactions]<-1
sampled_int_asym[lower.tri(sampled_int_asym,diag=TRUE)]<-t(sampled_int_asym)[lower.tri(sampled_int_asym,diag=TRUE)]
interaction_matrix<-sampled_int_asym

return(interaction_matrix)
}

### Creation of the matrices
#We create the adjacency matrices with different parameters. In particular, we vary the type of interactions (facilitation or competition), the percentage of interactions that we sample among all the ones that follow the criteria and the number of species. We plot the so-built adjacency matrix.

workdir="~/Phd/Elton John/Simul/V_1/data"

# Sp is 50 or 100
# We take niche_breath=20. This means that with 50 species, there are 345 possible interactions (each species can compete/facilitate with 10 species to the left and 10 to the right). With 100 species, there are 1449 species, and each species can compete/facilitate with 20 species to the left and 20 to the right.
# We take Nint = 10, 20, 50, 100

### Simulation 

sim_names=c(
  ### Faciliation
  "FacNint10Sp50",
  "FacNint20Sp50",
  "FacNint50Sp50",
  "FacNint100Sp50",

  "FacNint10Sp100",
  "FacNint20Sp100",
  "FacNint50Sp100",
  "FacNint100Sp100",

  ### Competition ###
  "CompNint10Sp50",
  "CompNint20Sp50",
  "CompNint50Sp50",
  "CompNint100Sp50",

  "CompNint10Sp100",
  "CompNint20Sp100",
  "CompNint50Sp100",
  "CompNint100Sp100"
)

nruns=length(sim_names)
set.seed(1123)


nruns=length(sim_names)
sim_params =list(
    Sp           = rep(rep(list(50,100),each=4),2),
    Nint         = rep(list(10,20,50,100),4),
    niche_breadth= rep(list(20),nruns),
    type         = rep(list("facilitation","competition"),each=8)
) 

adj_mats = sim_params %>% pmap(CreateInteractionMatrix)
adj_mats=  adj_mats %>% set_names(sim_names) 


lapply(sim_names, function(x){
  pdf(file=paste0(paste(workdir,"AdjMatrices/Plots",x,sep="/"),".pdf"))
  plot(adj_mats[[x]]==1,main=x,col=c('white', 'red'))
  dev.off()
  })

saveRDS(adj_mats[1:8],file=paste(workdir,"AdjMatrices/fac_adjmatrices.rds",sep="/"))
saveRDS(adj_mats[9:16],file=paste(workdir,"AdjMatrices/comp_adjmatrices.rds",sep="/"))


### Plot the Interaction networks


lapply(sim_names, function(x){
  plot(adj_mats[[x]]==1,main=x,col=c('white', 'red'))
  })


### Creating communities
#We create the communities following the parameters above. For each simulation, we take the interaction network created above.

outputdir=paste(workdir,"Communities",sep="/")

sim_names=c(
  ### Faciliation
  "FacNint10Sp50",
  "FacNint20Sp50",
  "FacNint50Sp50",
  "FacNint100Sp50",

  "FacNint10Sp100",
  "FacNint20Sp100",
  "FacNint50Sp100",
  "FacNint100Sp100",

  ### Competition ###
  "CompNint10Sp50",
  "CompNint20Sp50",
  "CompNint50Sp50",
  "CompNint100Sp50",

  "CompNint10Sp100",
  "CompNint20Sp100",
  "CompNint50Sp100",
  "CompNint100Sp100"
)

nruns=length(sim_names)
set.seed(1123)

lcomp=readRDS(file=paste(workdir,"AdjMatrices/comp_adjmatrices.rds",sep="/"))
lfac=readRDS(file=paste(workdir,"AdjMatrices/fac_adjmatrices.rds",sep="/"))


nruns=length(sim_names)

sim_params =list(
    niche_optima = rep(rep(list(seq(2, 98, length.out = 50),seq(2, 98, length.out = 100)),each=4),2),
    comp_inter   = c(rep(list(NA),8),lcomp),
    fac_inter    = c(lfac,rep(list(NA),8)),
    beta_comp    = rep(list(10),each=nruns),
    beta_fac     = rep(list(5),each=nruns),
    beta_env     = rep(list(1), each = nruns),
    beta_abun    = rep(list(0), each = nruns),
    K            = rep(list(40), each = nruns),
    competition  = rep(list("facilitation","competition"),each=8),
    intra_sp_com = rep(list(0), nruns),
    epochs       = rep(T,nruns)
) 
sim_data = sim_params %>% pmap(simulate_community) 
sim_data=sim_data %>% set_names(sim_names) 
saveRDS(sim_data,file=paste(outputdir,"sim_data.rds",sep="/"))

### Convergence validation
#Each boxplot represent the distribution of the shannon index over all communities at one time step. All the plots show convergence for any kind of parameter when we take 100 time steps.


lcomp=readRDS(file=paste(workdir,"AdjMatrices/comp_adjmatrices.rds",sep="/"))
lfac=readRDS(file=paste(workdir,"AdjMatrices/fac_adjmatrices.rds",sep="/"))
#sim_data=readRDS(file=paste(outputdir,"sim_data.rds",sep="/"))


alpha_div=function(tab){ ##tab with rows=epochs and columns=species
  a=diversity(tab,"shannon")
}

### nruns = nombre de simulations
### The following plots evolution of alpha-diversity (given by Shannon index applied to species counts) over epochs
### Each iteration corresponds to one simulation model

for (r in 1:nruns){
  hist=sim_data[[r]]$hist  ##contains the data for each epoch rows=epochs and columns=species
  df=do.call(rbind,lapply(hist,alpha_div)) #sites x years matrix
  pdf(file=paste0(paste(outputdir,"ConvPlots",sim_names[[r]],sep="/"),".pdf"))
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
  write.csv2(occur,file=paste(outputdir,paste0(x,".csv"),sep="/"))
})
