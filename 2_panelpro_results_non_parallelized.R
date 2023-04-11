#### Set Up ####
set.seed(777)
# load libraries
library(plyr) #need to load plyr before dplyr
library(dplyr)
library(PanelPRO)
library(tidyverse)

load("RObjects/familyinfo10Families.Rdata")


#### Full Families ####

fullCarrierRisk <- list()

#get proband carrier risk for full families
for(i in 1:length(families)){
  families[[i]]$MLH1 <- NA
  out <- PanelPRO::PanelPRO(families[[i]], genes="MLH1", cancers="Colorectal")
  id <- as.character(probandIDS[i])
  fullCarrierRisk[[i]] = out$posterior.prob[[id]]$estimate[2]
}

print("Ran PanelPRO on Full Families")


#### Full Families with Masked Information for Unaffected Individuals ####

# get outputs of PanelPRO for families with masked info
maskedCarrierRisk <- list()

#get proband carrier risk for full families with masked information about unaffected relatives
for(i in 1:length(families)){
  maskedFamilies[[i]]$MLH1 <- NA
  out <- PanelPRO::PanelPRO(maskedFamilies[[i]], genes="MLH1", cancers="Colorectal")
  id <- as.character(probandIDS[i])
  maskedCarrierRisk[[i]] = out$posterior.prob[[id]]$estimate[2]
}

print("Ran PanelPRO with unaffected family members info masked")


#### First Degree Families ####
firstDegreeCarrierRisk <- list()


#get proband carrier risk for first degree families nothing masked
for(i in 1:length(families)) {
  firstDegree[[i]]$MLH1 <- NA
  out <- PanelPRO::PanelPRO(firstDegree[[i]], genes="MLH1", cancers="Colorectal")
  id <- as.character(probandIDS[i])
  firstDegreeCarrierRisk[[i]] = out$posterior.prob[[id]]$estimate[2]
  
}


print("Ran PanelPRO on first degree families")
print(firstDegreeCarrierRisk)

#### Combine Results Into a Summary Table ####
summaryTableMLH1 <- cbind(format(unlist(fullCarrierRisk), scientific=FALSE), unlist(maskedCarrierRisk), unlist(firstDegreeCarrierRisk), unlist(probandMLH1Status))
summaryTableMLH1 <- as.data.frame(summaryTableMLH1)
summaryTable <- cbind(1:numberFamilies, summaryTableMLH1)
names(summaryTable) <- c("famID","fullCarrierRisk", "carrierRiskUnaffectedInfoMasked", "firstDegreeCarrierRisk", "probandMLH1Status")
print(head(summaryTable))
summaryTable$fullCarrierRisk <- as.numeric(summaryTable$fullCarrierRisk)
summaryTable$carrierRiskUnaffectedInfoMasked <- as.numeric(summaryTable$carrierRiskUnaffectedInfoMasked)
summaryTable$firstDegreeCarrierRisk <- as.numeric(summaryTable$firstDegreeCarrierRisk)


summaryTable <- cbind(summaryTable, unlist(numAffectedRelatives))
summaryTable <- cbind(summaryTable, unlist(numMLH1Relatives))
summaryTable <- cbind(summaryTable, unlist(numMLH1UnaffectedRelatives))

names(summaryTable) <- c("famID","fullCarrierRisk", "carrierRiskUnaffectedInfoMasked", "firstDegreeCarrierRisk", "probandMLH1Status", "numAffectedRelatives", "numMLH1Relatives","numMLH1UnaffectedRelatives")

View(summaryTable)


save(summaryTable, file = "RObjects/panelPROSummaryTable10FamiliesNoForEach.Rdata")

