#### Set Up ####

# load libraries
library(plyr) #need to load plyr before dplyr
library(dplyr)
library(PanelPRO)
library(tidyverse)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)

load("RObjects/familyinfo.Rdata")

cores <- 8
cl <- makeCluster(cores)
registerDoParallel(cl)

#### Full Families ####

outputsFull <- list ()
fullCarrierRisk <- list()

#get full outputs of PanelPRO for full families
outputsFull <- foreach(i = 1:length(families), .packages = c("PanelPRO")) %dopar% {
  out <- PanelPRO::PanelPRO(families[[i]], genes="MLH1", cancers="Colorectal")
  out
}

#get proband carrier risk for full families
fullCarrierRisk <- foreach(i = 1:length(families), .packages = c("PanelPRO")) %dopar% {
  out <- PanelPRO::PanelPRO(families[[i]], genes="MLH1", cancers="Colorectal")
  id <- as.character(probandIDS[i])
  out$posterior.prob[[id]]$estimate[2]
}

print("Ran PanelPRO on Full Families")


#### Full Families with Masked Information for Unaffected Individuals ####

# get outputs of PanelPRO for families with masked info
outputsMasked <- list ()
maskedCarrierRisk <- list()

#get full outputs of PanelPRO for full families with masked information about unaffected relatives
outputsMasked <- foreach(i = 1:length(families), .packages = c("PanelPRO")) %dopar% {
  out <- PanelPRO::PanelPRO(maskedFamilies[[i]], genes="MLH1", cancers="Colorectal")
  out
}

#get proband carrier risk for full families with masked information about unaffected relatives
maskedCarrierRisk <- foreach(i = 1:length(families), .packages = c("PanelPRO")) %dopar% {
  out <- PanelPRO::PanelPRO(maskedFamilies[[i]], genes="MLH1", cancers="Colorectal")
  id <- as.character(probandIDS[i])
  out$posterior.prob[[id]]$estimate[2]
}

#View(outputsMasked)
print("Ran PanelPRO with unaffected family members info masked")


#### First Degree Families ####
outputsFirstDegree <- list ()
firstDegreeCarrierRisk <- list()
#get full outputs of PanelPRO for full families with masked information about unaffected relatives
outputsFirstDegree <- foreach(i = 1:length(families), .packages = c("PanelPRO")) %dopar% {
  out <- PanelPRO::PanelPRO(firstDegree[[i]], genes="MLH1", cancers="Colorectal")
  out
}

#get proband carrier risk for full families with masked information about unaffected relatives
firstDegreeCarrierRisk <- foreach(i = 1:length(families), .packages = c("PanelPRO")) %dopar% {
  out <- PanelPRO::PanelPRO(firstDegree[[i]], genes="MLH1", cancers="Colorectal")
  id <- as.character(probandIDS[i])
  out$posterior.prob[[id]]$estimate[2]
}
print("Ran PanelPRO on first degree families")


#### Combine Results Into a Summary Table ####
summaryTableMLH1 <- cbind(format(unlist(fullCarrierRisk), scientific=FALSE), unlist(maskedCarrierRisk), unlist(firstDegreeCarrierRisk), unlist(probandMLH1Status))
summaryTableMLH1 <- as.data.frame(summaryTableMLH1)
summaryTable <- cbind(1:numberFamilies, summaryTableMLH1)
names(summaryTable) <- c("famID","fullCarrierRisk", "carrierRiskUnaffectedInfoMasked", "firstDegreeCarrierRisk", "probandMLH1Status")
print(head(summaryTable))
summaryTable$fullCarrierRisk <- as.numeric(summaryTable$fullCarrierRisk)
summaryTable$carrierRiskUnaffectedInfoMasked <- as.numeric(summaryTable$carrierRiskUnaffectedInfoMasked)
summaryTable$firstDegreeCarrierRisk <- as.numeric(summaryTable$firstDegreeCarrierRisk)


summaryTable <- cbind(summaryTable, numAffectedRelatives)
summaryTable <- cbind(summaryTable, numMLH1Relatives)
summaryTable <- cbind(summaryTable, numMLH1UnaffectedRelatives)
View(summaryTable)
stopCluster(cl)

save(summaryTable, file = "RObjects/panelPROSummaryTable.Rdata")

