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
library(foreach)
library(doParallel)



# Generate Families ----
numberFamilies <- 10000
families = list()
probandIDS = c()
probandMLH1Status = c()

cores <- 48
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


# registerDoParallel(cl)
# families <- foreach(i=1:numberFamilies, .packages="truncnorm") %dopar% {
#   # Cancers
#   cancers = "Colorectal"
#   # Genes
#   genes = "MLH1"
#   #family members
#   # Paternal aunts, paternal uncles
#   nSibsPatern =floor(rtruncnorm(n=2, mean=3, 3))
#   # Maternal aunts, maternal uncles
#   nSibsMatern = floor(rtruncnorm(n=2, mean=3, 3))
#   # Sisters and brothers
#   nSibs = floor(rtruncnorm(n=2, mean=3, 3))
#   # We make the assumption that the number of sons and daughters for the
#   # proband and all siblings, is the same. Nieces and nephews of the proband
#   # are not sampled separately.
#   nGrandchild = floor(rtruncnorm(n=2, mean=6, 2))
#   nChild = floor(rtruncnorm(n=2, mean=3, 2))
#   
#   # Simulate family using `PedUtils` code
#   fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nChild,
#                       PanelPRODatabase, genes, cancers,
#                       includeGeno = TRUE, includeBiomarkers = TRUE)
#   famDF = as.data.frame(fam)
#   for(j in 1:nrow(famDF)){
#     if(famDF[j,]$isAffAny == 0){
#       famDF[j,]$AgeAny = NA    
#     }
#     if(famDF[j,]$isAffCOL == 0){
#       famDF[j,]$AgeCOL = NA    
#     }
#   }
#   proband = famDF %>% filter(isProband==1)
#   probandIDS = c(probandIDS, proband$ID)
#   probandMLH1Status = c(probandMLH1Status, proband$MLH1)
#   # PanelPRO can be run on the simulated family
#   #out = PanelPRO:::PanelPRO11(fam)
#   families[[i]] = famDF
#   #outputs = c(outputs, out)
#   
# }


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

### ambry info both proband and relative ---- 
ambryInfo <- function(ped){
  proband <- ambryProbandInfo(ped)
  relatives <- ambryRelativeInfo(ped)
  family <- rbind.fill(proband, relatives)
  family <- affectedFamilyMembers(family)
  return(family)
}


## function to make a binary list of proband MLH1 carrier status ----
mlh1Probands <- function(ped){
  if(nrow(ped %>% filter(MLH1 ==1) %>% filter(isProband==1) > 0)){
    return(1)
  } 
  else{
    return(0)
  }
}




# Run Filters in parallel with mclapply ----
#cores <- 10
#cl <- makeCluster(cores)
clusterExport(cl, c("families"))
families <- mclapply(families, removeProbandStatus, mc.cores = cores, mc.preschedule=FALSE) #removes mlh1 status from probands
firstDegree <- mclapply(families, firstDegreeFamilyMembers, mc.cores = cores, mc.preschedule=FALSE) #filters to first degree
ambryFirstDegree <- mclapply(firstDegree, ambryInfo, mc.cores = cores, mc.preschedule = FALSE) #filters first degree info
mlh1ProbandsTotal <- mclapply(mlh1StatusFamilies, mlh1Probands, mc.cores = cores, mc.preschedule = FALSE) #determines which families have mlh1
stopCluster(cl)



outputsFull <- list ()
fullCarrierRisk <- list()
considered_fams = families
for(i in 1:numberFamilies){
  out = PanelPRO::PanelPRO(families[[i]], genes="MLH1", cancers="Colorectal")
  id = as.character(probandIDS[i])
  fullCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
  outputsFull[[i]] = out#DF
}


# #for first degree families
outputsFirstDegree <- list ()
firstDegreeCarrierRisk <- list()
considered_fams = firstDegree
for(i in 1:numberFamilies){
  out = PanelPRO::PanelPRO(firstDegree[[i]], genes="MLH1", cancers="Colorectal")
  id = as.character(probandIDS[i])
  firstDegreeCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
  outputsFirstDegree[[i]] = out#DF
}
outputsFirstDegreeAmbry <- list ()
firstDegreeAmbryCarrierRisk <- list()
#considers probands and their parents and affected family members
for(i in 1:numberFamilies){
  # browser()
  out = PanelPRO::PanelPRO(ambryFirstDegree[[i]], genes="MLH1", cancers="Colorectal", debug = TRUE)
  id = as.character(probandIDS[i])
  firstDegreeAmbryCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
  outputsFirstDegreeAmbry[[i]] = out#DF
}

affectedFamiliesLister <- function(f){
  affectedFamiliesList = list()
  for(i in 1:length(families)){
    fam = as.data.frame(f[[i]])
    fam = fam %>% filter(isAffAny==1) %>% filter(isProband==0)
    if(nrow(fam) > 0){
      affectedFamiliesList[[i]] = nrow(fam)
    }
    else{
      affectedFamiliesList[[i]] = 0}
  }
  return(affectedFamiliesList)
}

mlh1ProbandsTotal <- c()
probandMLH1 <- list()
for(i in 1:length(families)){
  if(nrow(mlh1StatusFamilies[[i]] %>% filter(MLH1 ==1) %>% filter(isProband==1) > 0)){
    mlh1ProbandsTotal <- c(mlh1ProbandsTotal, i)
    probandMLH1[[i]] = 1
  } 
  else{
    probandMLH1[[i]] = 0
  }}

# Generate a summary table ----
summaryTableMLH1 <- cbind(format(unlist(fullCarrierRisk), scientific=FALSE), unlist(firstDegreeCarrierRisk), unlist(firstDegreeAmbryCarrierRisk), unlist(affectedFamiliesLister(families)), unlist(affectedFamiliesLister(firstDegree)), unlist(probandMLH1))

summaryTableMLH1 <- as.data.frame(summaryTableMLH1)

names(summaryTableMLH1) <- c("Carrier Risk for Full Families", "Carrier Risk for First Degree Families", "Carrier Risk for Ambry Info", "Affected Individuals Full Familes", "Affected Individuals First Degree Families", "True MLH1 Carriers in Probands")

print(summaryTableMLH1)


summaryTable <- cbind(1:numberFamilies, summaryTableMLH1)
names(summaryTable) <- c("famNumber", "carrierRiskFullFamilies", "carrierRiskFirstDegreeFamilies", "carrierRiskAmbryInfo", "affectedIndividualsFullFamilies", "affectedIndividualsFirstDegreeFamilies", "probandsMLH1")
summaryTable

summaryTable$carrierRiskFullFamilies <- as.numeric(summaryTable$carrierRiskFullFamilies)
summaryTable$carrierRiskAmbryInfo <- as.numeric(summaryTable$carrierRiskAmbryInfo)
summaryTable$carrierRiskFirstDegreeFamilies <- as.numeric(summaryTable$carrierRiskFirstDegreeFamilies)
ggplot(data= summaryTable) +  geom_point(aes(x=famNumber, y=carrierRiskFullFamilies), color="Red") + 
  geom_point(aes(x=famNumber, y= carrierRiskFirstDegreeFamilies), color="Blue") + 
  geom_point(aes(x=famNumber, y=carrierRiskAmbryInfo), color="Green") +
  #geom_point(aes(x=famNumber, y=probandsMLH1), color = "Purple") +
  scale_y_log10() +
  labs(title= "MLH1 Carrier Risk", x= "Family ID", y="Carrier Risk")


ggplot(data= summaryTable) +  geom_point(aes(x=famNumber, y=carrierRiskFullFamilies, color="Carrier Risk for Full Families")) + 
  geom_point(aes(x=famNumber, y= carrierRiskFirstDegreeFamilies, color="Carrier Risk for First Degree Families")) + 
  geom_point(aes(x=famNumber, y=carrierRiskAmbryInfo, color="Carrier Risk Given Ambry Info")) +
  #geom_point(aes(x=famNumber, y=probandsMLH1, color = "Probands with MLH1")) +
  scale_y_log10() +
  labs(title= "MLH1 Carrier Risk", x= "Family ID", y="Carrier Risk (log 10 scale)")


# Median polish ----
fitMedianPolished <- rlm(carrierRiskFullFamilies ~ carrierRiskAmbryInfo, data = summaryTable)
summary(fitMedianPolished)

# finally, plot the data and the fitted curve
ggplot(data = summaryTable, aes(x= carrierRiskAmbryInfo, y= carrierRiskFullFamilies)) +
  geom_point() + 
  geom_line(aes(y = predict(fitMedianPolished)))+
  labs(x = "Ambry Carrier Risk First Degree Families", y= "Carrier Risk for Full Family Information", title = "Comparing Carrier Risk for Ambry and Full Family Information with Loess Regression")




# Loess ---- 
# next, use the loess() function to fit a model
fitLoess <- loess(carrierRiskFullFamilies ~ carrierRiskAmbryInfo, data = summaryTable, span = 0.2)
summary(fitLoess) #RSE is 0.00397 (seems good but I'm not sure)

# finally, plot the data and the fitted curve
ggplot(data = summaryTable, aes(x= carrierRiskAmbryInfo, y= carrierRiskFullFamilies)) +
  geom_point() + 
  geom_line(aes(y = predict(fitLoess)))+
  labs(x = "Ambry Carrier Risk First Degree Families", y= "Carrier Risk for Full Family Information", title = "Comparing Carrier Risk for Ambry and Full Family Information with Loess Regression") 




