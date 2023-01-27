set.seed(777)

# Load Libraries ----

setwd("../PedUtils/R")
files.sources = list.files()
sapply(files.sources, source) #loads in all the PedUtils functions
setwd("../../thesis-work/")
library(plyr) #need to load plyr before dplyr
library(dplyr)
library(truncnorm)
library(PanelPRO)
library(tidyverse)
library(MASS)
library(parallel)
library(parallel)


# Generate Families ----
numberFamilies <- 100
families = list()
probandIDS = c()
probandMLH1Status = c()

for (i in 1:numberFamilies) {
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


# Functions to Filter Family Info ----
## Functions to Remove MLH1 Status from Probands/Relatives ----
mlh1StatusFamilies <- families
removeProbandStatus <- function(ped){
  proband = ped %>% filter(isProband == 1)
  proband$MLH1 = NA
  relatives = ped %>% filter(isProband == 0)
  fam = rbind(proband, relatives)
  return(fam)
}
removeRelativeStatus <- function(ped){
  proband = ped %>% filter(isProband == 1)
  relatives = ped %>% filter(isProband == 0)
  relatives$MLH1 = NA
  fam = rbind(proband, relatives)
  return(fam)
}

## Function to Filter to First Degree Families ----
firstDegreeFamilyMembers <- function(ped){
  proband <- ped %>% filter(isProband==1)
  firstDegree = proband
  if(!is.null(proband)){
    for(i in  1:nrow(ped)){
      relative=ped[i,]
      #siblings
      if(relative$ID != proband$ID){
        if(!is.null(relative$MotherID) && !is.null(proband$MotherID) && !is.null(relative$FatherID) && !is.null(proband$FatherID) && (relative$MotherID==proband$MotherID) && (relative$FatherID==proband$FatherID)){
          firstDegree <- rbind(firstDegree, relative)
        }
        #parents
        else if(!is.null(proband$MotherID) && !is.null(proband$FatherID) && (relative$ID == proband$MotherID || relative$ID == proband$FatherID)){
          firstDegree <-rbind(firstDegree,relative)
        }
      }
    }
  }
  return(firstDegree)
}





## Functions to filter down to Ambry info ----
#### Proband info: sex, age, isDead==no, NA for riskMod and for interAge ###
ambryProbandInfo <- function(ped){
  proband = ped %>% filter(isProband==1)
  #proband$riskmod <- character(0)
  #proband$interAge <- character(0) #should these be character(0)? what does that mean?
  proband$isDead <- 0
  return(proband)
}

### relative info ###
#for simplicity, we will just say that sex can be inferred
#restricting to first degree family
#here genes is just MLH1
ambryRelativeInfo <- function(ped){
  relatives = ped %>% filter(isProband==0)
  #print(relatives)
  #relatives$CurAge <- NULL
  relatives$isDead <- NULL
  #relatives$riskMod <- character(0)
  #relatives$interAge <- character(0)
  #View(relatives)
  return(relatives)
}


## Function to filter to affected family members, proband, and parents ----
affectedFamilyMembers <- function(ped){
  proband <- ped %>% filter(isProband==1)
  mother <- ped %>% filter(ID == proband$MotherID)
  father <- ped %>% filter(ID == proband$FatherID)
  relatives <- ped %>% filter(isAffAny==1)
  relatives <- relatives %>% filter(isProband==0 && ID != proband$MotherID && ID != proband$FatherID) #to make sure we don't duplicate family members
  affectedFam <- rbind(proband, relatives, mother, father)
  return(affectedFam)
}




# Run Filters in parallel ----
cl <- makeCluster(24)
clusterExport(cl, c("families"))
families <- mclapply(families, removeProbandStatus, mc.cores = 24, mc.preschedule=FALSE)
firstDegree <- mclapply(families, firstDegreeFamilyMembers, mc.cores = 24, mc.preschedule=FALSE)
prob <- mclapply(firstDegree, ambryProbandInfo, mc.cores = cores, mc.preschedule = FALSE)
rels <- # NEED TO FINISH HERE LAUREN
  for(i in 1:length(firstDegree)){
    prob = ambryProbandInfo(firstDegree[[i]])
    rels = ambryRelativeInfo(firstDegree[[i]])
    famAm = rbind.fill(prob, rels)
    famAm = famAm %>% affectedFamilyMembers()
    ambryFirstDegree[[i]] = famAm
  }
stopCluster(cl)

