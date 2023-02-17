set.seed(777)

# Load Libraries ----

setwd("../PedUtils/R")
files.sources = list.files()
sapply(files.sources, source) #loads in all the PedUtils functions
setwd("../../thesis-work")
source("family_filters.R")
library(plyr) #need to load plyr before dplyr
library(dplyr)
library(truncnorm)
library(PanelPRO)
library(tidyverse)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)


generateFamilies <- function(mlh1Freq, numberFamilies){
  # Set MLH1 Frequency to 0.2 ----
  load("~/PanelPRO/data/PanelPRODatabase.rda")
  PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",1]] <- mlh1Freq
  PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",2]] <- mlh1Freq
  PanelPRODatabase$AlleleFrequency[["MLH1_anyPV",3]] <- mlh1Freq
  save(PanelPRODatabase, file = "~/PanelPRO/data/PanelPRODatabase.rda")
  
  
  # Generate Families ----
  families = list()
  probandIDS = c()
  probandMLH1Status = c()
  
  cores <- 25
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  families <- foreach(i=1:numberFamilies, .packages =c("truncnorm","tidyverse")) %dopar%{
    setwd("../PedUtils/R")
    files.sources = list.files()
    sapply(files.sources, source) #loads in all the PedUtils functions
    setwd("../../thesis-work/")
    # Cancers
    cancers = "Colorectal"
    # Genes
    genes = "MLH1"
    #family members
    # Paternal aunts, paternal uncles
    nSibsPatern =floor(rtruncnorm(n=2, mean=3, 3))
    # Maternal aunts, maternal uncles
    nSibsMatern = floor(rtruncnorm(n=2, mean=3, 3))
    # Sisters and brothers
    nSibs = floor(rtruncnorm(n=2, mean=3, 3))
    # We make the assumption that the number of sons and daughters for the
    # proband and all siblings, is the same. Nieces and nephews of the proband
    # are not sampled separately.
    nGrandchild = floor(rtruncnorm(n=2, mean=6, 2))
    nChild = floor(rtruncnorm(n=2, mean=3, 2))
    
    # Simulate family using `PedUtils` code
    fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nChild,
                        PanelPRODatabase, genes, cancers,
                        includeGeno = TRUE, includeBiomarkers = TRUE)
    famDF = as.data.frame(fam)
    for(j in 1:nrow(famDF)){
      if(famDF[j,]$isAffAny == 0){
        famDF[j,]$AgeAny = NA    
      }
      if(famDF[j,]$isAffCOL == 0){
        famDF[j,]$AgeCOL = NA    
      }
    }
    proband = famDF %>% filter(isProband==1)
    probandIDS = c(probandIDS, proband$ID)
    probandMLH1Status = c(probandMLH1Status, proband$MLH1)
    # PanelPRO can be run on the simulated family
    #out = PanelPRO:::PanelPRO11(fam)
    families[[i]] = famDF
    #outputs = c(outputs, out)
  }
  probandIDS <- foreach(i=1:numberFamilies, .packages="tidyverse") %dopar%{
    proband <- families[[i]] %>% filter(isProband==1)
    probandID <- proband$ID
  }
  probandMLH1Status <- foreach(i=1:numberFamilies, .packages = "tidyverse") %dopar%{
    proband <- families[[i]] %>% filter(isProband==1)
    probandMLH1Status <- proband$MLH1
  }
  mlh1StatusFamilies <- families
  families <- foreach(i=1:numberFamilies, .export='removeProbandStatus') %dopar%{
    families[[i]]=removeProbandStatus(families[[i]])
  }
  firstDegree <- foreach(i=1:numberFamilies, .export = 'firstDegreeFamilyMembers') %dopar%{
    firstDegreeFamilyMembers(families[[i]])
  }
  maskedFamilies <- foreach(i=1:numberFamilies, .export = 'maskedInfoAmbry') %dopar%{
    maskedInfoAmbry(families[[i]])
  }
  # clusterExport(cl, c("families")) #not sure if we need this if we generate families in the cluster
  # families <- mclapply(families, removeProbandStatus, mc.cores = cores, mc.preschedule=FALSE) #removes mlh1 status from probands
  # firstDegree <- mclapply(families, firstDegreeFamilyMembers, mc.cores = cores, mc.preschedule=FALSE) #filters to first degree
  # maskedFamilies <- mclapply(families, maskedInfoAmbry, mc.cores = cores, mc.preschedule = FALSE)
  # stopCluster(cl)
  
  print("Full Families from PedUtils with MLH1 Status for Probands")
  describeFamilies(mlh1StatusFamilies)
  
  print("Full Families from PedUtils with Unaffected Information Masked")
  describeFamilies(maskedFamilies)
  
  print("First Degree Families")
  describeFamilies(firstDegree)
  
  # cl <- makeCluster(cores)
  # registerDoParallel(cl)
  
  
  fullCarrierRisk <- foreach(i=1:numberFamilies, .packages="PanelPRO") %dopar% {
    out = PanelPRO::PanelPRO(families[[i]], genes="MLH1", cancers="Colorectal")
    id = as.character(probandIDS[i])
    out$posterior.prob[[id]]$estimate[2]
  }
  print("Ran PanelPRO for full carrier risks")
  
  # maskedFamiliesCarrierRisk <- foreach(i=1:numberFamilies, .packages="PanelPRO") %dopar% {
  #   out = PanelPRO::PanelPRO(maskedFamilies[[i]], genes="MLH1", cancers="Colorectal")
  #   id = as.character(probandIDS[i])
  #   out$posterior.prob[[id]]$estimate[2]
  # }
  # print("Ran PanelPRO with unaffected family members info masked")

  firstDegreeCarrierRisk <- foreach(i=1:numberFamilies, .packages="PanelPRO") %dopar% {
    out = PanelPRO::PanelPRO(firstDegree[[i]], genes="MLH1", cancers="Colorectal")
    id = as.character(probandIDS[i])
    out$posterior.prob[[id]]$estimate[2]
  }
  print("Ran PanelPRO on first degree families")
  stopCluster(cl)
 

  
}

families <- generateFamilies(0.2, 50)