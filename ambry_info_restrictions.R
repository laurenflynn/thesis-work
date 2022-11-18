set.seed(777)
setwd("/Users/flynnlauren7/Desktop/Internships Summer 2022/Bayes Mendel/PedUtils/R")
files.sources = list.files()
sapply(files.sources, source) #loads in all the PedUtils functions

#load libraries
library(dplyr)
library(truncnorm)
library(plyr)
library(PanelPRO)
library(ggplot2)


#generate families
families = list()
probandIDS = c()
#outputs = c()
for (i in 1:500) {
  # Cancers
  cancers = "Colorectal"
  # Genes
  genes = "MLH1"
  #family members
  # Paternal aunts, pat ernal uncles
  nSibsPatern =floor(rtruncnorm(n=2, mean=3, 2))
  # Maternal aunts, maternal uncles
  nSibsMatern = floor(rtruncnorm(n=2, mean=3, 2))
  # Sisters and brothers
  nSibs = floor(rtruncnorm(n=2, mean=3, 2))
  # We make the assumption that the number of sons and daughters for the
  # proband and all siblings, is the same. Nieces and nephews of the proband
  # are not sampled separately.
  nGrandchild = floor(rtruncnorm(n=2, mean=6, 2))
  nChild = floor(rtruncnorm(n=2, mean=3, 2))

  # Simulate family using `PedUtils` code
  fam = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nChild,
                      PanelPRODatabase, genes, cancers,
                      includeGeno = FALSE, includeBiomarkers = TRUE)
  famDF = as.data.frame(fam)
  proband = famDF %>% filter(isProband==1)
  probandIDS = c(probandIDS, proband$ID)
  # PanelPRO can be run on the simulated family
  #out = PanelPRO:::PanelPRO11(fam)
  families[[i]] = famDF
  #outputs = c(outputs, out)

}



## Function to filter down to Ambry info ##
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
  relatives$CurAge <- NULL
  relatives$isDead <- NULL
  #relatives$riskMod <- character(0)
  #relatives$interAge <- character(0)
  #View(relatives)
  return(relatives)
}


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


#returns all affected family members and proband (affected or otherwise)
affectedFamilyMembers <- function(ped){
  proband <- ped %>% filter(isProband==1)
  affectedFam = proband
  for(i in  1:nrow(ped)){
    relative=ped[i,]
    if( relative$ID != proband$ID && relative$isAffAny!=0){
      rbind(affectedFam, relative)
    }
  }
  return(affectedFam)
}


# first <- firstDegreeFamilyMembers(families[[3]])
# # View(first)
# 
# ambryInfoFam <- ambryRelativeInfo(first)
# ambryInfoProband <- ambryProbandInfo(first)
# # View(ambryInfoFam)
# ambryInfo <- rbind.fill(ambryInfoFam, ambryInfoProband)
# # View(ambryInfo)
# affected <- affectedFamilyMembers(ambryInfo)
# View(affected)

#filter down to first degree families
firstDegreeFamilies = list()
for(i in 1:length(families)){
  first <- as.data.frame(firstDegreeFamilyMembers(families[[i]]))
  firstDegreeFamilies[[i]] = first
}


#restrict to Ambry Info and affected family members only
ambryInfoAffectedFamilies = list()
for(i in 1:length(families)){
  ambryInfoFam <- ambryRelativeInfo(firstDegreeFamilies[[i]])
  ambryInfoProband <- ambryProbandInfo(firstDegreeFamilies[[i]])
  ambryInfo <- rbind.fill(ambryInfoFam, ambryInfoProband)
  affected <- affectedFamilyMembers(ambryInfo)
  #if just the proband is left, take away the mother id and father id
  proband <- affected %>% filter(isProband == TRUE)
  if(!(proband$MotherID %in% affected$ID)){
    proband$MotherID <- NA
  }
  if(!(proband$FatherID %in% affected$ID)){
    proband$FatherID <- NA
  }
  ambryInfoAffectedFamilies[[i]] <- affected
  
}


#### testing how many families have affected relatives
affectedFamilies = 0 
affectedFAMILIES = list()
j = 1
for(i in 1:length(families)){
  f = firstDegreeFamilies[[i]]
  f = f %>% filter(isAffAny==1) %>% filter(isProband==0)
  if(nrow(f)>0){
    affectedFamilies = affectedFamilies + 1
    affectedFAMILIES[[j]] = f
    j = j + 1
  }
}
affectedFamilies #25 families have affected family members, 0 probands


#generating outputs for first degree families with full information
outputsFull <- list ()
fullCarrierRisk <- list()
for(i in 1:length(families)){
  out = PanelPRO::PanelPRO(firstDegreeFamilies[[i]], , genes="MLH1", cancers="Colorectal")
  id = as.character(probandIDS[i])
  fullCarrierRisk[[i]] <- out$posterior.prob[[id]]$estimate[2]
  outputsFull[[i]] = out#DF
}

#generating outputs for first degree families with Ambry information
outputsAmbry <- list ()
ambryRisks <- list()
for(i in 1:length(families)){
  out = PanelPRO::PanelPRO(ambryInfoAffectedFamilies[[i]], genes="MLH1", cancers="Colorectal")
  id = as.character(probandIDS[i])
  outRisk = out$posterior.prob[[id]]$estimate[2]
  ambryRisks[[i]] = outRisk
  #outDF = as.data.frame(out)
  outputsAmbry[[i]] = out#DF
}


mlh1CarrierEstimatesFirstDegree <- outputsFull[[i]][["posterior.prob"]][["21"]]


plot(fullCarrierRisk, ambryRisks)
plot(ambryRisks, fullCarrierRisk)



fullRisk <- unlist(fullCarrierRisk)
ambryRisk <- unlist(ambryRisks)
riskDF <- as.data.frame(cbind(ambryRisk,fullRisk))
riskDF2 <- as.data.frame(cbind(ambryRisk,fullRisk, ambryInfoAffectedFamilies, firstDegreeFamilies))
ggplot(data=riskDF) + geom_point(aes(x = ambryRisk, y = fullRisk)) + labs(x="Risks from Ambry Restricted Info", y="Risk from First Degree Families", title="Mapping Carrier Risk from Ambry Restricted Info to First Degree Family Information")
                                                                      
removed_outlier <- riskDF %>% filter(fullRisk < 0.2)

outlier = riskDF2 %>% filter(fullRisk > 0.2)
ggplot(data=removed_outlier) + geom_point(aes(x = ambryRisk, y = fullRisk)) + labs(x="Risks from Ambry Restricted Info", y="Risk from First Degree Families", title="Mapping Carrier Risk from Ambry Restricted Info to First Degree Family Information")

#including outlier--linear regression
summary(lm(riskDF$fullRisk~riskDF$ambryRisk))

#excluding outlier--linear regression
summary(lm(removed_outlier$fullRisk~removed_outlier$ambryRisk))














